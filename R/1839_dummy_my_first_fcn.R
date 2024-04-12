#install package dependencies
if(!require('devtools')){
   install.packages('devtools')
}

#Check that dependencies have been installed
devtools::dev_package_deps(dependencies = T)

#This is a function created for the admiral pharmaverse package
# dummy example to print when run
my_first_fcn <-  function(x){
  print("Welcome to the admiral family!")
}

my_first_fcn()
