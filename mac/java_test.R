# Diagnostics for testing you have the correct version of Java on your computer. 
# Tests are sequential (fromhttps://stackoverflow.com/questions/14915898/rnetlogo-function-nlstart-fails-to-launch-gui).     

# test 1  

.jinit()
.jnew( "java/awt/Point", 10L, 10L )
# f <- .jnew("java/awt/Frame","Hello")
# .jcall(f,,"setVisible",TRUE)
# # cat("\n","\n",geterrmessage(),"\n","\n")
# # stop("Failed Test 1",call.=F)
# 
# 
# # test 2
# component <- .jnull()
# component <- .jcast(component, new.class = "java/awt/Component")
# message <- .jnew("java/lang/String","This is a JOptionPane test from rJava.")
# message <- .jcast(message, new.class = "java/lang/Object")
# title <- .jnew("java/lang/String","Test")
# type <- .jnew("java/lang/Integer", as.integer(2))
# f <- .jnew("javax/swing/JOptionPane")
# .jcall(f,"showMessageDialog", component, message, title, .jsimplify(type))
# cat("\n","\n",geterrmessage(),"\n","\n")
# stop("Failed Test 2",call.=F)
# 
# 
# # test 3
# .jcall("java/lang/System", "S", "getProperty", "java.vm.version")
# .jcall("java/lang/System", "S", "getProperty", "java.vm.name")
# .jcall("java/lang/System", "S", "getProperty", "java.vm.info")
# .jcall("java/lang/System", "S", "getProperty", "java.runtime.version")
# .jcall("java/lang/System", "S", "getProperty", "sun.arch.data.model")

# test 4
.jcall("java/lang/System", "S", "getProperty", "java.awt.headless")
Sys.getenv("NOAWT")



