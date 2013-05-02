thinslice_doc
=============

This model simulates the ignition-extinction hysteresis of the catalytic combustion takes place in a thin slice Diesel Oxidation Catalyst (DOC).
A Diesel Oxidation Catalyst is a pollution control device that treats the diesel exhaust gas before it is discharged into the atmosphere.
A DOC is a catalyst system with Pt as the catalyst deposited on alumina washcoat (wc). The catalyst support is cordierite monolith structure.

There are two sets of models:
- pde_1wc_thin_hc_y.m pairs with wc_thin_hc_y.m     ---- these are the 1D gas phase - 2D solid phase model
- pde_1wc_thin_hc_y2D.m pairs with wc_thin_hc_y2D.m ---- these are the 1D gas phase - 2D solid phase model

The models consider pollutants such as CO, hydrocarbon (assume propane), NOx which are the typical pollutants found in the untreated diesel exhaust gas.
The models consider a thin slice (5mm) monolith structure. This is to approximate the structure as a differential reactor and allow detailed study.
The models use the oxygen compression mechanisms as proposed by Salomons et al (2007) to account for the high density oxygen phase during the extinction process.

Generally, the:
- pde_1wc_thin_hc m-files comprise the PDEs that describe the dispersion-advection in the monolith channel, and the diffusion-reaction in the porous washcoat
- wc_thin_hc m-files comprise parameters, conversion from PDEs to ODEs and generation of plots.

The m-files will generate plots of ignition-extinction hysteresis.

Shoot me a message if you need any clarification on the code.
