### Copyright (C) 1998, 1999
### Justin Lokhorst <jlokhors@stats.adelaide.edu.au>
### Berwin A. Turlach <bturlach@stats.adelaide.edu.au>
### Bill Venables <wvenable@stats.adelaide.edu.au>
### $Id: First.lib.template,v 1.18 1999/12/08 23:35:43 bturlach Exp $
### --> ../COPYRIGHT for more details

.First.lib <- function(lib, pkg)
{
    library.dynam("lasso2",pkg,lib)

    cat("R Package to solve regression problems while imposing",
        "\t an L1 constraint on the parameters. Based on S-plus Release 2.1",
        "Copyright (C) 1998, 1999",
        "Justin Lokhorst   <jlokhors@stats.adelaide.edu.au>",
        "Berwin A. Turlach <bturlach@stats.adelaide.edu.au>",
        "Bill Venables     <wvenable@stats.adelaide.edu.au>\n",

        "Copyright (C) 2002",
        "Martin Maechler <maechler@stat.math.ethz.ch>\n", sep="\n")
}

