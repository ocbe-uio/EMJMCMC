simplify.formula<-function(fmla,names)
{
fmla.proc<-as.character(fmla)[2:3]
fobserved <- fmla.proc[1]
fmla.proc[2]<-stringi::stri_replace_all(str = fmla.proc[2],fixed = " ",replacement = "")
fmla.proc[2]<-stringi::stri_replace_all(str = fmla.proc[2],fixed = "\n",replacement = "")
fparam <- names[which(names %in% stringi::stri_split_fixed(str = fmla.proc[2],pattern = "+",omit_empty = F)[[1]] )]
return(list(fparam = fparam,fobserved = fobserved))
}
