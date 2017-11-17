#A function for when you load your script into someone else's computer, so it checks if present, install if not and load them all

uptodate<-function (list_of_packages){
  to_install<-list_of_packages[!list_of_packages %in% installed.packages()]
  if (length(to_install>0)){
  lapply(to_install,install.packages(), character.only=T)
  }
  lapply(list_of_packages,require, character.only=T)
}
