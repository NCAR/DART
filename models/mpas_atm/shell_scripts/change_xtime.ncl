; change xtime in AmB.nc

infn="AmB.nc"
a=addfile(infn,"w")
t=a->xtime
print("change_xtime.ncl")

ts=systemfunc("grep config_start_time namelist.atmosphere | awk '{print $3}'")
ts = str_strip(str_sub_str(ts,"'",""))
;print(ts+"")
tnew=tochar(ts)
t(0,:18) = (/tnew/)
a->xtime = (/t/)

;system("ls -l "+infn)
system("ncdump -v xtime "+infn+" | tail -2")
print("")
