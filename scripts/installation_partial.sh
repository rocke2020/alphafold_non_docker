### **Download chemical properties to the common folder**
export alphafold_path="$(pwd)/alphafold-2.3.1"
# wget -q -P $alphafold_path/alphafold/common/ https://git.scicore.unibas.ch/schwede/openstructure/-/raw/7102c63615b64735c4941278d92b554ec94415f8/modules/mol/alg/src/stereo_chemical_props.txt


### **Apply OpenMM patch**
cd ~/anaconda3/envs/af2/lib/python3.8/site-packages/ && patch -p0 < $alphafold_path/docker/openmm.patch