#!/bin/sh

export R_HOME="/Library/Frameworks/R.framework/Resources"
export R_ARCH=""
export R_LIBS="/Library/Frameworks/R.framework/Versions/3.2/Resources/library"
export R_LIBS_USER="~/Library/R/3.2/library"
'/Library/Java/JavaVirtualMachines/1.6.0.jdk/Contents/Home/bin/java' -cp '/Library/Frameworks/R.framework/Versions/3.2/Resources/library/rJava/java/boot' -Drjava.class.path='/Library/Frameworks/R.framework/Versions/3.2/Resources/library/rJava/jri/JRI.jar:/Library/Frameworks/R.framework/Versions/3.2/Resources/library/iplots/java/iplots.jar:/Library/Frameworks/R.framework/Versions/3.2/Resources/library/JGR/java/JGR.jar:/Library/Frameworks/R.framework/Resources/etc/classes:/Library/Frameworks/R.framework/Resources/etc/classes/classes.jar' -Drjava.path='/Library/Frameworks/R.framework/Versions/3.2/Resources/library/rJava' -Dmain.class=org.rosuda.JGR.JGR  -Djgr.load.pkgs=yes  -Dr.arch= -Xdock:icon='/Library/Frameworks/R.framework/Versions/3.2/Resources/library/JGR/icons/JGR.icns' -Dcom.apple.mrj.application.apple.menu.about.name='JGR' -Xdock:name='JGR' -Dapple.laf.useScreenMenuBar=true -Dcom.apple.macos.useScreenMenuBar=true -Xrs -Xss5m RJavaClassLoader 

