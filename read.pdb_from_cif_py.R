require(reticulate)
read.pdb_from_cif_py <- function(cif_path=NULL,
                        use_youngest_alt=TRUE,
                        python_path="/Users/sakuma/PycharmProjects/mmcif/venv/bin/python3.8"){
  use_python(python = python_path)
  pdbecif=import("pdbecif")
  reader=pdbecif$mmcif_io$CifFileReader()
  cif_out=reader$read(file_path = cif_path,output='cif_dictionary')
  
  #
  cl=match.call()
  
  #make pdb object 
  pdb=list()
  #class(pdb)="pdb"
  
  # type eleno elety  alt resid chain resno insert     x       y      z o     b
  # segid elesy charge
  n_atom = length(cif_out[[1]]$`_atom_site`$group_PDB)
  atom=data.frame(type=as.character(cif_out[[1]]$`_atom_site`$group_PDB),
                  eleno=as.numeric(cif_out[[1]]$`_atom_site`$id),
                  elety=as.character(cif_out[[1]]$`_atom_site`$label_atom_id),
                  alt=as.character(cif_out[[1]]$`_atom_site`$label_alt_id),
                  resid=as.character(cif_out[[1]]$`_atom_site`$label_comp_id),
                  chain=as.character(cif_out[[1]]$`_atom_site`$label_asym_id),
                  resno=as.numeric(cif_out[[1]]$`_atom_site`$label_seq_id),
                  insert=as.character(cif_out[[1]]$`_atom_site`$pdbx_PDB_ins_code),
                  x=as.numeric(cif_out[[1]]$`_atom_site`$Cartn_x),
                  y=as.numeric(cif_out[[1]]$`_atom_site`$Cartn_y),
                  z=as.numeric(cif_out[[1]]$`_atom_site`$Cartn_z),
                  o=as.numeric(cif_out[[1]]$`_atom_site`$occupancy),
                  b=as.numeric(cif_out[[1]]$`_atom_site`$B_iso_or_equiv),
                  segid=rep(NA,n_atom),
                  elesy=as.character(cif_out[[1]]$`_atom_site`$type_symbol),
                  charge=as.character(cif_out[[1]]$`_atom_site`$pdbx_formal_charge),
                  stringsAsFactors = F)
  
  pdb$atom = atom
  
  xyz=vector(length = 3*n_atom)
  ind_x = seq(1,3*n_atom-2,3)
  ind_y = seq(2,3*n_atom-1,3)
  ind_z = seq(3,3*n_atom,3)
  xyz[ind_x]=as.numeric(cif_out[[1]]$`_atom_site`$Cartn_x)
  xyz[ind_y]=as.numeric(cif_out[[1]]$`_atom_site`$Cartn_y)
  xyz[ind_z]=as.numeric(cif_out[[1]]$`_atom_site`$Cartn_z)
  pdb$xyz=xyz
  
  pdb$atom$alt[which(pdb$atom$alt==".")]=NA
  pdb$atom$insert[which(pdb$atom$insert=="?")]=NA
  pdb$atom$resno[which(is.na(pdb$atom$resno))]=0
  
  class(pdb)=c("pdb","sse")
  
  pdb$helix = NULL
  pdb$sheet = NULL
  pdb$seqres = NULL
  pdb$calpha = NULL
  pdb$remark = NULL
  pdb$call = cl
  
  
  
  
  
  pdb_new=pdb
  if (use_youngest_alt){
    resnos=unique(pdb$atom$resno[which(!is.na(pdb$atom$alt))])
    if(length(resnos)>0){
      indexes=c()
      for (resno in resnos){
        ind_resno=which(pdb$atom$resno==resno)
        youngest_alt=na.omit(unique(pdb$atom$alt[ind_resno]))[1]
        print(paste0("Using the youngest alt=",youngest_alt," for resno=",resno))
        ind_resno_youngest_alt=unique(which( ( is.na(pdb$atom$alt)| pdb$atom$alt==youngest_alt )& pdb$atom$resno==resno))
        indexes=c(indexes,ind_resno_youngest_alt)
      }
      ind_noalt=which(is.na(pdb$atom$alt) & !is.element(pdb$atom$resno,resnos))
      index_all=sort(c(indexes,ind_noalt))
      pdb_new$atom=pdb$atom[index_all,]
      pdb_new$xyz=pdb$xyz[atom2xyz(index_all)]
    }
  }
  return(pdb_new)
  
}