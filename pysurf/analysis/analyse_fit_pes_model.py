from pysurf.workflow import engine

workflow = engine.create_workflow(
    "analyse_fit_pes_model",
    """
spp = spp_analysis("spp.inp")
spp_calc = spp_calc("spp.inp", 5, 3, ['energy'])
sampler = sampler("samplerinp")
crds = crds_from_sampler(sampler, 200)
energies = get_energies(spp, crds)
energies_calc = get_energies(spp_calc, crds)
array2D_scale(energies, [0, 1, 2], 27.2114)
array2D_scale(energies_calc, [0, 1, 2], 27.2114)
x = array_extract_column(crds, 0)
data = array2D_add_column(energies, x, append=False)
data_calc = array2D_add_column(energies_calc, x, append=False)
datasorted = array2D_sort(data)
datasorted_calc = array2D_sort(data_calc)
plot = setup_lineplot("plot.inp")
color = standard_colors()
dashed = linestyle_dashed()
style_dashed = combine_plotstyles(color, dashed)
add_plot(plot, datasorted_calc)
add_plot(plot, datasorted, style=style_dashed)
show_plot(plot)
""",
)

workflow.run()
