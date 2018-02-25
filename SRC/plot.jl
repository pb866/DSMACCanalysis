### Output file definition and initialisation ###

println("plot data...")
# Define default output file name and directory
ofile = ARGS[2]
if ofile == ""  ofile = join(label,"_")  end
if dirname(ofile) == ""
  ofile = normpath(joinpath(Base.source_dir(),"../../save/results/",ofile))
end
if lowercase(splitext(ofile)[end]) != ".pdf"  ofile *= ".pdf"  end
# Initilise multipage output pdf
pdffile = pdf.PdfPages(ofile)
# Define night-time shading
nights = def_night(rates,pltnight[2])

#––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––#

### Plotting loop ###

# Loop over different plot types
for n = 1:length(icase)
  if what[n]=="fluxes"
    # Loop over species and scenarios for flux plots
    for spc in plotdata[n]  for case in icase[n]
      if sfig=="on"
        sfile = join([what[n],spc,label[case]],"_")*".pdf"
        sfile = joinpath(Base.source_dir(),"../FIG/"*replace(sfile,"/","-"))
      else
        sfile = ""
      end
      # Load plot data from ropa analysis
      src, snk, src_rev, snk_rev =
        load_plotdata(spc,case,sources,sinks,concs,unit[n],llim=lims[1],ulim=lims[2])
      # Define time format
      modtime = specs[case][Symbol(t_frmt)]
      # Output flux plots
      fig = plot_data(spc,label[case],modtime,src,snk,unit[n],
            nights[case],pltnight,t_frmt,sfile)
      if fig != nothing  pdffile[:savefig](fig)  end
      # if fig != nothing  fig[:show]()  end
      # Output revised flux plots, if major fluxes have been removed
      fig = plot_data(spc,label[case],modtime,src_rev,snk_rev,unit[n],
            nights[case],pltnight,t_frmt,sfile)
      if fig != nothing  pdffile[:savefig](fig)  end
      # if fig != nothing  fig[:show]()  end
      # input("Next picture?")
    end  end
  elseif what[n]=="stack"
    colstyle = ["source","sink"]
    # Plot stack plots of species concentrations for all cases
    for spc in plotdata[n]  for i = 1:length(icase[n])
      if sfig=="on"
        if length(spc) > 3
          spc_list = join(spc[1:3],"+")*"+more"
        else
          spc_list = join(spc[1:3],"+")
        end
        sfile = join([what[n],spc_list,label[icase[n][i]]],"_")*".pdf"
        sfile = joinpath(Base.source_dir(),"../FIG/"*replace(sfile,"/","-"))
      else
        sfile = ""
      end
      ylines, ystack = get_stackdata(spc,icase[n][i],specs,unit[n])
      lt = make_plots.sel_ls(cs=colstyle[i], nc=1:length(ylines), nt=icase[n][i])
      # Define time format
      modtime = specs[icase[n][i]][Symbol(t_frmt)]
      fig = make_plots.plot_stack(modtime,ylines,ystack,label[icase[n][i]],spc,
                                  unit[n],lt,nights[icase[n][i]],pltnight,t_frmt,sfile)
      pdffile[:savefig](fig)
    end  end
  elseif what[n] == "specs"
    # Define night-time shading for current plot section
    curr_night = pltnight
    night = []
    try night = nights[icase[n][1]]
      # Test night-times are the same, if not, omit shading for this section
      for i = 2:length(icase[n])
        if nights[1] != nights[i]
          curr_night[2] = 0.0
          println("Warning! Different night-times in plot section $n.")
          println("Night-time shading switched of for this case.")
          break
        end
      end
    catch; curr_night = ["w", 0.0]
    end
    # Define time format
    modtime = specs[icase[n][1]][Symbol(t_frmt)]
    # Plot line plots of species concentrations for all cases
    for case in plotdata[n]
      if sfig=="on"
        sfile = join([what[n],join(case,"+"),join(label,"+")],"_")*".pdf"
        sfile = joinpath(Base.source_dir(),"../FIG/"*replace(sfile,"/","-"))
      else
        sfile = ""
      end
      fig = make_plots.lineplot(modtime,specs,label,what[n],unit[n],
                                icase[n],case,night,curr_night,t_frmt,sfile)
      pdffile[:savefig](fig)
    end
  elseif what[n] == "rates"
    # Define night-time shading for current plot section
    curr_night = pltnight
    night = []
    try night = nights[icase[n][1]]
      # Test night-times are the same, if not, omit shading for this section
      for i = 2:length(icase[n])
        if nights[1] != nights[i]
          curr_night[2] = 0.0
          println("Warning! Different night-times in plot section $n.")
          println("Night-time shading switched of for this case.")
          break
        end
      end
    catch; curr_night = ["w", 0.0]
    end
    # Define time format
    modtime = rates[icase[n][1]][Symbol(t_frmt)]
    # Plot line plots of reaction rates for all cases
    for case in plotdata[n]
      if sfig=="on"
        sfile = join([what[n],join(case,"+"),join(label,"+")],"_")*".pdf"
        sfile = joinpath(Base.source_dir(),"../FIG/"*replace(sfile,"/","-"))
      else
        sfile = ""
      end
      fig = make_plots.lineplot(modtime,rates,label,what[n],unit[n],
                                icase[n],case,night,curr_night,t_frmt,sfile)
      pdffile[:savefig](fig)
    end
  end
end
# Close multipage pdf file
pdffile[:close]()
