require(reticulate)
read.cif_py <- function(cif_path=NULL,
                                 python_path="/Users/sakuma/PycharmProjects/mmcif/venv/bin/python3.8"){
  use_python(python = python_path)
  pdbecif=import("pdbecif")
  reader=pdbecif$mmcif_io$CifFileReader()
  cif_out=reader$read(file_path = cif_path,output='cif_dictionary')
  class(cif_out)="cif_dict"
  return(cif_out)
}