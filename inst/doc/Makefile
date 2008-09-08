all:
	if [ -f Manual-wo-help.pdf ]; then 	\
	   texi2dvi --pdf Manual-vignette.tex;	\
	   rm *.log *.aux;			\
	   rm -rf auto;				\
	   rm .build.timestamp;			\
	fi
