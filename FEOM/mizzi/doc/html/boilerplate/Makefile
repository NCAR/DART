
include ../../../mkmf/mkmf.template

ALL = html_namelist html_function html_subroutine html_type html_functoc

default: $(ALL)

html_namelist : html_namelist.f90
	$(FC) $(FFLAGS) $< -o $@

html_function: html_function.f90
	$(FC) $(FFLAGS) $< -o $@

html_subroutine: html_subroutine.f90
	$(FC) $(FFLAGS) $< -o $@

html_type: html_type.f90
	$(FC) $(FFLAGS) $< -o $@

html_functoc: html_functoc.f90
	$(FC) $(FFLAGS) $< -o $@


clean:
	rm -f $(ALL) fort.22
