program ps_rand_local

use random_seq_mod, only : random_seq_type, init_random_seq, random_uniform

implicit none

integer :: num_sets, level, obs_kind, num, num_done
real :: err_var, bot_lat, top_lat, bot_lon, top_lon, lat, lon
type(random_seq_type) :: r

! Initialize the random sequence
call init_random_seq(r)

! Set up constants
num_sets = 1
level = -1
obs_kind = 3

! Open an output file and write header info
open(unit = 20, file = 'ps_rand.out')
write(20, *) 'set_def.out'
write(20, *) num_sets

write(*, *) 'input the number of observations'
read(*, *) num
write(20, *) num

write(*, *) 'input the obs error variance'
read(*, *) err_var

write(*, *) 'input a lower bound on latitude -90 to 90'
read(*, *) bot_lat
write(*, *) 'input an upper bound on latitude -90 to 90'
read(*, *) top_lat
write(*, *) 'input a lower bound on longitude '
read(*, *) bot_lon
write(*, *) 'input an upper bound on longitude '
read(*, *) top_lon

! Simple error check to let people know limits
if(top_lat <= bot_lat .or. top_lon <= bot_lon) then
   write(*, *) 'lat lon range error'
   stop
endif

num_done = 0
do while(num_done < num)
! Find a lat/lon pair in the box
   lat = bot_lat + random_uniform(r) * (top_lat - bot_lat)
   lon = bot_lon + random_uniform(r) * (top_lon - bot_lon)
   write(*, *) lat, lon

! Found one, output it
   num_done = num_done + 1
   write(20, *) err_var
   write(20, *) -1
   write(20, *) level
   write(20, *) lon
   write(20, *) lat
   write(20, *) obs_kind

end do

end program ps_rand_local
