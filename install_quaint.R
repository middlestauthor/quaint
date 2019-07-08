packages_needed <- c("devtools","roxygen2")

for (i in 1:length(packages_needed)){
  if(!(packages_needed[i] %in% installed.packages())){install.packages(packages_needed[i])}
}

for (i in 1:length(packages_needed)){
  library( packages_needed[i], character.only = TRUE)
}


document()
setwd('../')
install('quaint')


