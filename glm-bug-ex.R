##- From: Bjoern Reineking <bjorn@oesa.ufz.de>
##- To: maechler@stat.math.ethz.ch
##- Subject: predict.gl1ce(,type="response") yields values outside [0,1]
##- Date: Tue, 04 Jun 2002 15:49:40 +0100

##- Hallo Martin,

##- zunaechst noch einmal herzlichen Dank fuer Deine spontane Hilfe gestern!
##- Ich habe nun noch einmal die Lasso2-Bibliothek ausprobiert, und das
##- Ergebnis ist immer noch merkwuerdig (siehe beigefuegte Grafik):
##- a) die vorhergesagten Werte auf der "response"-Skala sind z.T. kleiner als 0.
##- b) die vorhergesagten Werte auf "response"-Skala zeigen keinen klaren
##- funktionalen Zusammenhang mit denen auf der "link"-Skala
##- Bei der "normalen" logistischen Regression sind beide Phaenomene nicht
##- vorhanden.
##-
##- Herzliche Gruesse,
##- Bjoern


## Als Code habe ich verwendet:

library(lasso2)
data(esoph)
modEso <- formula(cbind(ncases, ncontrols) ~ agegp + tobgp * alcgp)
glm.E  <- glm(modEso, data = esoph, family = binomial())
gl1c.E <- gl1ce(modEso, data = esoph, family = binomial())
png("gl1ce.png", width=480, height=960)
par(mfrow=c(2,1))
plot(predict(glm.E, type="link"), predict(glm.E, type="response"))
plot(predict(gl1c.E, type="link"), predict(gl1c.E, type="response"))
dev.off()

system("xv gl1ce.png &")# same plot as Björn ..

##-  Ich arbeite unter Windows NT:
version
##- 	  _
##-  platform i386-pc-mingw32
##-  arch     i386
##-  os       mingw32
##-  system   i386, mingw32
##-  status
##-  major    1
##-  minor    5.0
##-  year     2002
##-  month    04
##-  day      29
##-  language R


## DESCRIPTION-File meiner Lasso2-Version:

##  Package: lasso2
##  Version: 1.0-2
##  Author: Justin Lokhorst, Bill Venables and
##    Berwin Turlach <berwin@maths.uwa.edu.au>; first port to R: Martin Maechler
##  Maintainer: Martin Maechler <maechler@stat.math.ethz.ch>
##  Title: L1 constrained estimation aka `lasso'
##  Description: Routines and documentation for solving regression problems
##    while imposing an L1 constraint on the estimates, based on
##    the algorithm of Osborne et al. (1998)
##  Depends:
##  License: GPL version 2 or later
##  URL: http://www.maths.uwa.edu.au/~berwin/software/lasso.html
##  Built: R 1.4.1;  Win32;  2002-04-12 11:20:56


## Beispiel Grafik
## 	  [Content-Type: application/octet-stream; name="gl1ce.png"]
##
## --> ~/R/Pkgs/lasso2/glm-bug-gl1ce.png
##     ---------------------------------
