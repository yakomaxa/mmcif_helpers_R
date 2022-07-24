require(reticulate)
require(stringr)
read.cif_py <- function(cif_path=NULL,
                                 python_path="",
                                 ignore_atom_site=T,remove_underbar=T,verbose=F
                                 ){
  use_python(python = python_path)
  pdbecif=import("pdbecif")
  reader=pdbecif$mmcif_io$CifFileReader()
  if (ignore_atom_site){
    cif_out=reader$read(file_path = cif_path,output='cif_dictionary',ignore="_atom_site")
  }else{
    cif_out=reader$read(file_path = cif_path,output='cif_dictionary')
  }
  
  if(remove_underbar){
    n=names(cif_out[[1]])
    nn=str_remove(n,"_")
    names(cif_out[[1]])=nn
    if (verbose){
      cat("read.cif_py: The underbar at the head of first layers' element names are removed\n")
    }
  }
  
  class(cif_out)="cif_dict"
  return(cif_out)
}