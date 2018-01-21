# Get Julian time from model output
jtime=output["specs"][1][:TIME]
# Calculate correct time for first entry in array (set to 0 by DSMACC)
jtime[1]=2jtime[2]-jtime[3]
# Get Julian day and time of the day in seconds (jsec)
jday=Int.(jtime÷3600÷24)
jsec=jtime.-jday.⋅3600.⋅24
# Get latitude in radian from model output
lat = deg2rad.(output["specs"][1][:LAT])

# Calculate hourly angles in radians
w = deg2rad.(15.0*(jsec./3600-12.))

# Calculate declination
φ = 2.0π.⋅(jday-1)/365
dec = (0.006918 - 0.399912 .⋅ cos.(φ) + 0.070257 .⋅ sin.(φ)
      - 0.006758 .⋅ cos.(2.0 .⋅ φ) + 0.000907 .⋅ sin.(2.0 .⋅ φ)
      - 0.002697 .⋅ cos.(3.0 .⋅ φ) + 0.00148 .⋅ sin.(3.0 .⋅ φ))

# Calculate sza
χ = acos.(sin.(dec) .⋅ sin.(lat) .+ cos.(dec) .⋅ cos.(lat) .⋅ cos.(w))
sza = rad2deg.(χ)


sunrise = find([sza[i]>90 && sza[i+1]<90 for i = 1:length(sza)-1])
sunset  = find([sza[i]<90 && sza[i+1]>90 for i = 1:length(sza)-1])

if sunrise[1] < sunset[1]  insert!(sunset,1,1)  end
if sunrise[end] < sunset[end]  push!(sunset,1)  end
if length(sunrise) != length(sunset)
  println("Unequal number of sunrises and sunsets. Script stopped."); exit()
end
using PyPlot

fig, ax = subplots()
plot(output["time"],output["rates"][1][Symbol("NO2-->O+NO")])
for n = 1:length(sunrise)
  ax[:axvspan](output["time"][sunset[n]], output["time"][sunrise[n]],
    alpha=0.2, color="black", lw=0)
end
plt[:show]()
