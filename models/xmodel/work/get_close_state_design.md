# `get_close_state` Design for xmodel

## Overview

`get_close_state` uses a two-phase approach: a single monolithic spatial search over the
full combined state vector (all models and all lags), followed by per-model distance
modification. This avoids maintaining per-model `get_close_type` objects and keeps lag
handling implicit.

---

## Phase 1 — Monolithic spatial search

```fortran
call loc_get_close_obs(gc, base_loc, base_type, locs, qtys, num_close, close_ind, distances)
```

`gc` is built by filter over the full concatenated `locs` array. All models and all lag
copies are included, so the spatial search returns correct candidates for the entire state
vector in a single call.

---

## Phase 2 — Per-model distance modification

For each base model, extract the subset of `close_ind` that belongs to that model,
call its distance modifier, then write the (possibly changed) distances back:

```fortran
do i = 1, num_base_models
    call extract_model_subset(i, num_close, close_ind, distances, &
                              sub_count, sub_close_ind, sub_distances, sub_local_indx)
    call <model>_modify_close_state_distances(base_loc, base_type, locs, qtys, &
                                              sub_local_indx, sub_count, sub_close_ind, &
                                              sub_distances, state_handle)
    call write_back_distances(i, sub_count, sub_close_ind, distances)
enddo
```

Identifying which `close_ind` entries belong to model `i` is cheap:

```fortran
local_pos = close_ind(j) - model_offsets(i) + 1
if (local_pos >= 1 .and. local_pos <= model_sizes(i)) ! belongs to model i
```

No new data structures are needed.

---

## The `modify_close_state_distances` interface

A new routine added to each `model_mod.f90`:

```fortran
subroutine modify_close_state_distances(base_loc, base_type, locs, qtys, &
                                        local_indx, num_close, close_ind, &
                                        distances, state_handle)
   type(location_type), intent(in)    :: base_loc
   integer,             intent(in)    :: base_type
   type(location_type), intent(in)    :: locs(:)
   integer,             intent(in)    :: qtys(:)
   integer(i8),         intent(in)    :: local_indx(:)   ! indices local to this model
   integer,             intent(in)    :: num_close
   integer,             intent(in)    :: close_ind(:)    ! indices into locs
   real(r8),            intent(inout) :: distances(:)    ! modified in place
   type(ensemble_type), optional, intent(in) :: state_handle
```

A model that needs no special behaviour provides a one-line no-op default. Models that
need to suppress cross-level or cross-domain increments set specific `distances` entries
to a large value.

`modify_close_state_distances` is added to the `ROUTINES` array in
`generate_assim_model_mod.sh` so it is renamed and imported like all other model_mod
routines.

---

## Lag "distance" scaling

For a smoother, lag copies of a model are earlier time windows. Should older lags have difference localization?. More lag => greater "distance" from the obs? 

### Lag models in the generated code

The generator emits a dedicated modifier for each lag model. It first delegates to the
base model's modifier, then applies a lag-number penalty:

Note the modifier is the get_close modifier mentioned above. Thinking about multiple models, multiple lags.

```fortran
subroutine camfvlag2_modify_close_state_distances(...)
   ! Apply base model logic first
   call camfv_modify_close_state_distances(...)
   ! Inflate distances for lag 2
   distances(1:num_close) = distances(1:num_close) * lag_distance_factor(2)
end subroutine
```

### `lag_distance_factor` — runtime namelist

A module-level array read from the xmodel namelist:

```fortran
real(r8) :: lag_distance_scale = 1.0_r8   ! namelist parameter
real(r8), allocatable :: lag_distance_factor(:)   ! computed at init

! At init, for lag_num = 1..NLAGS:
lag_distance_factor(lag_num) = lag_distance_scale ** lag_num
```

With `lag_distance_scale = 1.5`:

| Lag | Factor |
|-----|--------|
| 1   | 1.5    |
| 2   | 2.25   |
| 3   | 3.375  |
| 4   | 5.06   |

Setting `lag_distance_scale = 1.0` disables lag scaling (all lags treated equally).

### Physical interpretation

Distance inflation narrows the effective localisation window in observation space.
A state variable at lag 4 is only updated by very nearby observations. This approximates
the expected decay of cross-time covariance with lag without requiring explicit covariance
modelling.

---

## Summary of advantages over per-model `get_close_type`

| Property | Per-model `gc` | Distance modifier |
|---|---|---|
| Spatial search called once | no | yes |
| Separate `get_close_type` per model built at init | yes | no |
| Lag handling | explicit index replication | free — lags already in full `gc` |
| Interface models implement | full `get_close_state` | simple `modify_close_state_distances` |
| Models need global offsets | no | no |
| Lag decay supported | manual | via `lag_distance_scale` namelist |

**Note** I think models do need global offsets because x,y,z from dart_index is the global dart index.

---

## Implementation steps

1. Add `modify_close_state_distances` to `ROUTINES` in `generate_assim_model_mod.sh`
2. Add a default no-op `modify_close_state_distances` to a shared utility (used by models that don't override it)
3. Add `lag_distance_scale` to the xmodel namelist and emit `lag_distance_factor` initialisation in `static_init_assim_model`
4. Replace the `get_close_state` stub in the generator with the two-phase implementation
5. Emit lag-specific modifier subroutines for each lag model in the generator
