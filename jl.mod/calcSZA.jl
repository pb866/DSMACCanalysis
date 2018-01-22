# Get Julian day and time of the day in seconds (jsec)
jday = Dates.dayofyear.(specs[1][:JTIME])
# jd = Dates.value.(Dates.Second.(Dates.Day.(jday)))
jsec = Dates.value.(Dates.Second.(Dates.Hour.(specs[1][:JTIME]))) .+
       Dates.value.(Dates.Second.(Dates.Minute.(specs[1][:JTIME]))) .+
       Dates.value.(Dates.Second.(specs[1][:JTIME]))

# Get latitude in radian from model output
lat = deg2rad.(specs[1][:LAT])

# Calculate hourly angles in radians
w = deg2rad.(15.0*(jsec./3600-12.))

# Calculate declination
φ = 2.0π.⋅(jday-1)/365
dec = (0.006918 - 0.399912cos.(φ) + 0.070257sin.(φ)
      - 0.006758cos.(2.0φ) + 0.000907sin.(2.0φ)
      - 0.002697cos.(3.0φ) + 0.00148sin.(3.0φ))

# Calculate sza
χ = acos.(sin.(dec) .⋅ sin.(lat) .+ cos.(dec) .⋅ cos.(lat) .⋅ cos.(w))
sza = rad2deg.(χ)


sunrise = find([sza[i]>90 && sza[i+1]<90 for i = 1:length(sza)-1])
sunset  = find([sza[i]<90 && sza[i+1]>90 for i = 1:length(sza)-1])

if sunrise[1] < sunset[1]  unshift!(sunset,1)  end
if sunrise[end] < sunset[end]  push!(sunrise,length(sza))  end
if length(sunrise) != length(sunset)
  println("Unequal number of sunrises and sunsets. Script stopped."); exit()
end
using PyPlot

fig, ax = subplots()
plot(rates[1][:TIME],rates[1][Symbol("NO2-->O+NO")])
for n = 1:length(sunrise)
  ax[:axvspan](rates[1][:TIME][sunset[n]], rates[1][:TIME][sunrise[n]],
    alpha=0.2, color="black", lw=0)
end
plt[:show]()
