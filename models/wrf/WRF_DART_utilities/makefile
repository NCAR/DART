#
#  pgi compilers
#
#.F.o:
#	pgf90 -g -Mfree -I/usr/local/netcdf-3.5.0-pgi/include -c $*.F
#dart_tf_wrf.o : wrf_data_module.o
#dart_tf_wrf: dart_tf_wrf.o 
#	pgf90 -g -o dart_tf_wrf dart_tf_wrf.o /usr/local/netcdf-3.5.0-pgi/lib/libnetcdf.a
#wrf_meta_data: wrf_meta_data.o 
#	pgf90 -g -o wrf_meta_data wrf_meta_data.o /usr/local/netcdf-3.5.0-pgi/lib/libnetcdf.a


#
# intel compilers
#
dart_tf_wrf.o : wrf_data_module.o dart_tf_wrf.f90
	ifc -FR -c dart_tf_wrf.f90 -I../netcdf/include
dart_tf_wrf: dart_tf_wrf.o 
	ifc -o dart_tf_wrf dart_tf_wrf.o ../netcdf/lib/libnetcdf.a
wrf_meta_data: wrf_meta_data.o
	ifc -o wrf_meta_data wrf_meta_data.o ../netcdf/lib/libnetcdf.a
wrf_meta_data.o: wrf_meta_data.f90
	ifc -c wrf_meta_data.f90 -I../netcdf/include
wrf_data_module.o: wrf_data_module.f90
	ifc -FR -c wrf_data_module.f90 -I../netcdf/include

