cmake_minimum_required ( VERSION 2.8 )
project ( talk_fosdem_2013 )


find_package ( LATEX )

if ( LATEX_COMPILER )

   include ( UseLATEX.cmake )

   add_latex_document ( talk.tex
	IMAGES 
               ambrosys.png
               granular_system.png
               logo.pdf
               lorenz.jpg
               neuron.jpg
               solar_system.jpg
               lorenz_perf.pdf
               lyap.jpg
               turbulence.jpg
               lattices.jpg
               lorenz.jpg
	DEFAULT_PDF
	MANGLE_TARGET_NAMES
	)

endif ( LATEX_COMPILER )
