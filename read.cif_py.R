require(reticulate)
read.cif_py <- function(cif_path=NULL,
                                 python_path="",
                                 ignore_atom_site=T
                                 ){
  use_python(python = python_path)
  pdbecif=import("pdbecif")
  reader=pdbecif$mmcif_io$CifFileReader()
  if (ignore_atom_site){
    cif_out=reader$read(file_path = cif_path,output='cif_dictionary',ignore="_atom_site")
  }else{
    cif_out=reader$read(file_path = cif_path,output='cif_dictionary')
  }
  class(cif_out)="cif_dict"
  return(cif_out)
}