# Diagnostics for testing you have the correct version of Java on your computer.
# Tests are sequential (from https://stackoverflow.com/questions/14915898/rnetlogo-function-nlstart-fails-to-launch-gui).

# test 1
install.packages("rJava", repos = "https://cran.r-project.org/", type="source"); library(rJava)
.jinit()
.jnew( "java/awt/Point", 10L, 10L )
f <- .jnew("java/awt/Frame","Hello")
.jcall(f,,"setVisible",TRUE)
t1err <- geterrmessage()

# # test 2
component <- .jnull()
component <- .jcast(component, new.class = "java/awt/Component")
message <- .jnew("java/lang/String","This is a JOptionPane test from rJava.")
message <- .jcast(message, new.class = "java/lang/Object")
title <- .jnew("java/lang/String","Test")
type <- .jnew("java/lang/Integer", as.integer(2))
f <- .jnew("javax/swing/JOptionPane")
.jcall(f,"showMessageDialog", component, message, title, .jsimplify(type))
t2err <- geterrmessage()

# # test 3
.jcall("java/lang/System", "S", "getProperty", "java.vm.version")
.jcall("java/lang/System", "S", "getProperty", "java.vm.name")
.jcall("java/lang/System", "S", "getProperty", "java.vm.info")
.jcall("java/lang/System", "S", "getProperty", "java.runtime.version")
.jcall("java/lang/System", "S", "getProperty", "sun.arch.data.model")
t3err <- geterrmessage()

# test 4
.jcall("java/lang/System", "S", "getProperty", "java.awt.headless")
Sys.getenv("NOAWT")
t4err <- geterrmessage()

errorlist <- function(){
if(geterrmessage()==t1err){stop("Failed Test 1 — Headless exception \n \n Wrong awt GUI support for Java/rJava",call.=T)}
if(geterrmessage()==t2err){stop("Failed Test 2 — Invalid method name for RcallMethod (unable to open dialog box)",call.=T)}
if(geterrmessage()==t3err){stop("Failed Test 3 — Old version of Java. \n \n Download latest version \n \n > https://www.oracle.com/technetwork/java/javase/downloads/index-jsp-138363.html.",call.=T)}
if(geterrmessage()==t4err){stop("Failed Test 4 — ",call.=T)}
}



