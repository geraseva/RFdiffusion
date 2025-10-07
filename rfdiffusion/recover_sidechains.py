import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np

import pickle


from rfdiffusion.chemical import num2aa, aa2num, aa_123, seq2chars
from rfdiffusion.util import rigid_from_3_points

def nsplit(*x):
    return [i.split() for i in x]
mass = {'H': 1, 'C': 12, 'N': 14, 'O': 16, 'S': 32, 'P': 31, 'M': 0}
bb = "N CA C O H H1 H2 H3 O1 O2"
unknown_mapping=nsplit(bb + " CB")
mapping = {
        "ALA":  nsplit(bb + " CB"),
        "CYS":  nsplit(bb, "CB SG"),
        "SEC":  nsplit(bb, "CB SEG"),
        "ASP":  nsplit(bb, "CB CG OD1 OD2"),
        "GLU":  nsplit(bb, "CB CG CD OE1 OE2"),
        "ASX":  nsplit(bb, "CB CG ND1 ND2 OD1 OD2 HD11 HD12 HD21 HD22"),
        "GLX":  nsplit(bb, "CB CG CD OE1 OE2 NE1 NE2 HE11 HE12 HE21 HE22"),
        "PHE":  nsplit(bb, "CB CG CD1 HD1", "CD2 HD2 CE2 HE2", "CE1 HE1 CZ HZ"),
        "GLY":  nsplit(bb),
        "HIS":  nsplit(bb, "CB CG", "CD2 HD2 NE2 HE2", "ND1 HD1 CE1 HE1"),
        "HIH":  nsplit(bb, "CB CG", "CD2 HD2 NE2 HE2", "ND1 HD1 CE1 HE1"),     # Charged Histidine.
        "ILE":  nsplit(bb, "CB CG1 CG2 CD CD1"),
        "LYS":  nsplit(bb, "CB CG CD", "CE NZ HZ1 HZ2 HZ3"),
        "LEU":  nsplit(bb, "CB CG CD1 CD2"),
        "MET":  nsplit(bb, "CB CG SD CE"),
        "ASN":  nsplit(bb, "CB CG ND1 ND2 OD1 OD2 HD11 HD12 HD21 HD22"),
        "PRO":  nsplit(bb, "CB CG CD"),
        "HYP":  nsplit(bb, "CB CG CD OD"),
        "GLN":  nsplit(bb, "CB CG CD OE1 OE2 NE1 NE2 HE11 HE12 HE21 HE22"),
        "ARG":  nsplit(bb, "CB CG CD", "NE HE CZ NH1 NH2 HH11 HH12 HH21 HH22"),
        "SER":  nsplit(bb, "CB OG HG"),
        "THR":  nsplit(bb, "CB OG1 HG1 CG2"),
        "VAL":  nsplit(bb, "CB CG1 CG2"),
        "TRP":  nsplit(bb, "CB CG CD2", "CD1 HD1 NE1 HE1 CE2", "CE3 HE3 CZ3 HZ3", "CZ2 HZ2 CH2 HH2"),
        "TYR":  nsplit(bb, "CB CG CD1 HD1", "CD2 HD2 CE2 HE2", "CE1 HE1 CZ OH HH"),
        "DA": nsplit("N9 C4",     
                     "C8 N7 C5",
                     "P OP1 OP2 O5' H5T O3' H3T O1P O2P",
                     "C5' O4' C4'",
                     "C3' C2' C1'",
                     "C2 N3",
                     "C6 N6 N1"),
        "DG": nsplit("N9 C4",
                     "C8 N7 C5",
                     "P OP1 OP2 O5' H5T O3' H3T O1P O2P",
                     "C5' O4' C4'",
                     "C3' C2' C1'",
                     "C2 N2 N3",
                     "C6 O6 N1"),
        "DC": nsplit("N1 C6",
                     "P OP1 OP2 O5' H5T O3' H3T O1P O2P",
                     "C5' O4' C4'",
                     "C3' C2' C1'",
                     "N3 C2 O2",
                     "C5 C4 N4"),
        "DT": nsplit("N1 C6",
                     "P OP1 OP2 O5' H5T O3' H3T O1P O2P",
                     "C5' O4' C4'",
                     "C3' C2' C1'",   
                     "N3 C2 O2",
                     "C5 C4 O4 C7 C5M"),
        "A":  nsplit("N9 C4",
                     "C8 N7 C5",
                     "P OP1 OP2 O5' H5T O3' H3T O1P O2P",
                     "C5' O4' C4'",
                     "C3' C2' O2' C1'", 
                     "C2 N3",
                     "C6 N6 N1"),
        "G":  nsplit("N9 C4",
                     "C8 N7 C5",
                     "P OP1 OP2 O5' H5T O3' H3T O1P O2P",
                     "C5' O4' C4'",
                     "C3' C2' O2' C1'", 
                     "C2 N2 N3",
                     "C6 O6 N1"),
        "C":  nsplit("N1 C6",
                     "P OP1 OP2 O5' H5T O3' H3T O1P O2P",
                     "C5' O4' C4'",
                     "C3' C2' O2' C1'",
                     "N3 C2 O2",
                     "C5 C4 N4"),
        "U":  nsplit("N1 C6",
                     "P OP1 OP2 O5' H5T O3' H3T O1P O2P",
                     "C5' O4' C4'",
                     "C3' C2' O2' C1'",
                     "N3 C2 O2",
                     "C5 C4 O4 C7 C5M")} 
unknown_types=['P4']
pseudoatom_types = {
        "ALA":  ['P4'],
        "CYS":  ['P5','C5'],
        "SEC":  ['P5','C5'],
        "ASP":  ['P5','Qa'],
        "GLU":  ['P5','Qa'],
        "ASX":  ['P5','P4'],
        "GLX":  ['P5','P4'], 
        "PHE":  ['P5','SC4','SC4','SC4'],
        "GLY":  ['P5'],
        "HIS":  ['P5','SC4','SP1','SP1'],
        "ILE":  ['P5','AC1'],
        "LYS":  ['P5','C3','P1'],
        "LEU":  ['P5','AC1'],
        "MET":  ['P5','C5'],
        "ASN":  ['P5','P5'],
        "PRO":  ['P5','AC2'],
        "GLN":  ['P5','P4'],
        "ARG":  ['P5','N0','Qd'],
        "SER":  ['P5','P1'],
        "THR":  ['P5','P1'],
        "VAL":  ['P5','AC2'],
        "TRP":  ['P5','SC4','SP1','SC4','SC4'],
        "TYR":  ['P5','SC4','SC4','SP1'],
        "DA": ["TN0", "TNa", "Q0", "SN0","SC2", "TA2", "TA3"],
        "DG": ["TN0", "TNa", "Q0","SN0","SC2", "TG2", "TG3"],
        "DC": ["TN0", "Q0","SN0","SC2", "TY2", "TY3"],
        "DT": ["TN0", "Q0","SN0","SC2", "TT2", "TT3"],
        "A": ["TN0", "TNa", "Q0","SN0","SNda", "TA2", "TA3"],
        "G": ["TN0", "TNa", "Q0","SN0","SNda", "TG2", "TG3"],
        "C": ["TN0", "Q0","SN0","SNda", "TY2", "TY3"],
        "U": ["TN0", "Q0","SN0","SNda", "TT2", "TT3"],
} # from martini_v2.1

pseudoatom_radii = {
    'P5': 4.7, 'AC1': 4.7,'C5': 4.7, 'SP1': 4.7, 
    'N0': 4.7, 'AC2': 4.7,'C3': 4.7, 'P1': 4.7,
    'Qa': 4.7, 'P4': 4.7,'Qd': 4.7, 'SC4': 4.7,
    "Q0": 4.7, "SN0": 4.3, "SC2": 4.3, "SNda": 4.3, 
    "TN0": 3.2, "TA2": 3.2, "TA3": 3.2, "TG2": 3.2, "TG3": 3.2, 
    "TY2": 3.2, "TY3": 3.2, "TT2": 3.2, "TT3": 3.2, "TNa": 3.2    
} # from martini_v2.1

pseudoatom_weights = {
    'P5': 72.0, 'AC1': 72.0, 'C5': 72.0, 'SP1': 45.0, 
    'N0': 72.0, 'AC2': 72.0, 'C3': 72.0, 'P1': 72.0,
    'Qa': 72.0, 'P4': 72.0, 'Qd': 72.0, 'SC4': 45.0, 
    "Q0": 72.0, "SN0": 45.0, "SC2": 45.0, "SNda": 45.0, 
    "TN0": 45.0, "TA2": 45.0, "TA3": 45.0, "TG2": 45.0, "TG3": 45.0, 
    "TY2": 45.0, "TY3": 45.0, "TT2": 45.0, "TT3": 45.0, "TNa": 45.0
} # from martini_v2.1

with open('/home/domain/geraseva/geraseva_data/auxiliary_potential/masif_martini/datasets/ideal_coords.pkl', 'rb') as f:
    ideal_coords = pickle.load(f)

def martinize(seq, atoms, coords):
    ps_coords=[]
    ps_types=[]
    for i, (aa, ass, xyzs) in enumerate(zip(seq, atoms, coords)):

        av=[]
        for a, xyz in zip(ass, xyzs):
            for j, m in enumerate(mapping.get(aa,unknown_mapping)):
                if len(av)<=j:
                    av.append([[0,0,0],0])
                if a in m:
                    av[j][0][0]+=xyz[0]*mass[a[0]]
                    av[j][0][1]+=xyz[1]*mass[a[0]]
                    av[j][0][2]+=xyz[2]*mass[a[0]]
                    av[j][1]+=mass[a[0]]
        av=[[ps[0][0]/ps[1],ps[0][1]/ps[1],ps[0][2]/ps[1]] for ps in av if ps[1]>0]
        
        ps_coords.append(av)
        ps_types.append(pseudoatom_types.get(aa, unknown_types)[:len(av)])
    return ps_coords, ps_types

init_N = torch.tensor([-0.5272, 1.3593, 0.000]).float()
init_CA = torch.zeros_like(init_N)
init_C = torch.tensor([1.5233, 0.000, 0.000]).float()
norm_N = init_N / (torch.norm(init_N, dim=-1, keepdim=True) + 1e-5)
norm_C = init_C / (torch.norm(init_C, dim=-1, keepdim=True) + 1e-5)
cos_ideal_NCAC = torch.sum(norm_N*norm_C, dim=-1) # cosine of ideal N-CA-C bond angle


def get_O_from_3_points(xyz, bond=1.24, eps=1e-8):
   
    N=xyz[:,0,:]
    CA=xyz[:,1,:]
    C=xyz[:,2,:]

    N=torch.cat((N[1:,:],N[:1,:]), dim=0)

    v=((C-CA)/(torch.norm(C-CA, dim=-1, keepdim=True)+eps) +
       (C-N)/(torch.norm(C-N, dim=-1, keepdim=True)+eps))

    O=C+v/(torch.norm(v, dim=-1, keepdim=True)+eps)*bond

    xyz_out=torch.zeros_like(xyz)
    xyz_out[:,:3,:]=xyz[:,:3,:]
    xyz_out[:,3,:]=O

    return xyz_out

def probability(preds, labels, reduction='mean'):
    preds=torch.sigmoid(preds)
    if labels.sum()<=0:
        preds=1-preds
    if reduction=='none':
        return preds
    elif reduction=='mean':
        return preds.mean()
    elif reduction=='sum':
        return preds.sum()


class GetMartiniSidechains:

    def __init__(self, binderlen=-1, seq_model_type='protein_mpnn'):

        submodule_path='/'.join(__file__.split('/')[:-3])
        import sys
        sys.path.append(submodule_path)

        import LigandMPNN
        from LigandMPNN.model_utils import ProteinMPNN
        from LigandMPNN.data_utils import restype_str_to_int, restype_1to3, alphabet, featurize, element_list
        
        restype_3to1={restype_1to3[x]: x for x in restype_1to3.keys()}
        self.renumber_aa_mpnn2rf=torch.tensor([restype_str_to_int[restype_3to1.get(x,'X')] for x in num2aa], dtype=int)
        self.renumber_aa_rf2mpnn=torch.tensor([aa2num[aa_123.get(x,'UNK')] for x in alphabet], dtype=int)
        self.element_dict = dict(zip(element_list, range(1, len(element_list))))

        path_to_LigandMPNN=LigandMPNN.__path__._path[0]

        self.binderlen=binderlen

        self.device='cuda' if torch.cuda.is_available() else 'cpu'

        # Initialize LigandMPNN model
        if seq_model_type == "protein_mpnn":
            checkpoint_path = f"{path_to_LigandMPNN}/model_params/proteinmpnn_v_48_020.pt"
        elif seq_model_type == "ligand_mpnn":
            checkpoint_path = f"{path_to_LigandMPNN}/model_params/ligandmpnn_v_32_010_25.pt"
        elif seq_model_type == "soluble_mpnn":
            checkpoint_path = f"{path_to_LigandMPNN}/model_params/solublempnn_v_48_020.pt"

        seq_checkpoint = torch.load(checkpoint_path, map_location=self.device)
        
        self.seq_model = ProteinMPNN(
            node_features=128,
            edge_features=128,
            hidden_dim=128,
            num_encoder_layers=3,
            num_decoder_layers=3,
            k_neighbors=seq_checkpoint["num_edges"],
            device=self.device,
            atom_context_num=seq_checkpoint.get('atom_context_num',1),
            model_type=seq_model_type,
            ligand_mpnn_use_side_chain_context=0,
        )
        self.seq_model.load_state_dict(seq_checkpoint["model_state_dict"])
        self.seq_model.to(self.device)
        self.seq_model.eval()
        self.featurize=featurize
        print('Load LigandMPNN model')

        self.recover_sc=None

    def init_recover_sc(self):
        # Initialize BB2Martini
        
        #self.recover_sc=BB2Martini(chains= ['_p1'] if self.binderlen<0 else ['_p1','_p2'])
        self.recover_sc=BB2MartiniModule(chains= ['_p1'] if self.binderlen<0 else ['_p1','_p2'],
                                         device=self.device)
        self.recover_sc.eval()

    def run_LigandMPNN(self, xyz, seq, seq_mask, ligand_xyz=None, ligand_aatypes=None):

        L=xyz.shape[0]

        xyz=get_O_from_3_points(xyz)

        input_dict = {}
        
        input_dict["X"] = xyz[:,:4,:] # L*4*3 (bb atoms)  ? normalize
        input_dict["mask"] = torch.ones([ L]).to(self.device)

        if seq==None:
            input_dict["S"] = torch.full(( L),20,dtype=int).to(self.device) # encoded sequence
            input_dict["chain_mask"] = torch.full(( L),True,dtype=bool).to(self.device)
            raise AttributeError
        else:
            input_dict["S"]=seq[:,self.renumber_aa_rf2mpnn].argmax(-1).squeeze().detach()
            input_dict["chain_mask"] = ~seq_mask.squeeze().detach()

        if ligand_xyz==None:
            input_dict["Y"] = torch.zeros([1, 3]).to(self.device)
            input_dict["Y_t"] = torch.zeros([1]).to(self.device)
            input_dict["Y_m"] = torch.zeros([1]).to(self.device)
        else:
            input_dict["Y"] = ligand_xyz.to(self.device).detach()
            input_dict["Y_t"] = torch.tensor([self.element_dict.get(x,0) for x in ligand_aatypes], 
                                             dtype=torch.int32, device=self.device)
            input_dict["Y_m"] = torch.ones_like(input_dict["Y_t"])

        input_dict["R_idx"] = torch.arange(L).to(self.device) # L resnums

        if self.binderlen>0:
            input_dict["chain_labels"] = torch.cat((torch.zeros((self.binderlen)),
                                              torch.ones((L-self.binderlen))),0).to(self.device)  # L Chain indices
        else:
            input_dict["chain_labels"]=torch.zeros((L)).to(self.device)

        feature_dict = self.featurize(input_dict,
                                      number_of_ligand_atoms=(self.seq_model.features.atom_context_num 
                                                              if self.seq_model.model_type=='ligand_mpnn' 
                                                              else 1) ,
                                      model_type=self.seq_model.model_type)

        feature_dict["batch_size"]=1
        feature_dict["temperature"] = 0.1
        feature_dict["bias"] = torch.zeros((1,L,21)).to(self.device)
        feature_dict["randn"]=torch.randn((1,L)).to(self.device)
        feature_dict["symmetry_residues"] = [[]]
        feature_dict["symmetry_weights"]=[[]]

        output_dict = self.seq_model.score(feature_dict, use_sequence=True)

        return output_dict


    def get_aa_probs(self, xyz, seq, seq_mask, ligand_xyz=None, ligand_aatypes=None):

        output_dict=self.run_LigandMPNN(xyz, seq, seq_mask, ligand_xyz, ligand_aatypes)
        probs=torch.nn.functional.softmax(output_dict['logits'], dim=-1)
        probs=probs[0,:,self.renumber_aa_mpnn2rf]
        probs[seq_mask.squeeze()]=seq[seq_mask.squeeze()]
        print(f'ProteinMPNN prediction: { seq2chars(torch.argmax(probs, dim=-1).tolist())}')
        return probs.contiguous()

    def bb2martini(self, xyz, seq):
     
        feature_dict = {}

        if self.binderlen<0:
            feature_dict['sequence_p1']=seq
            feature_dict['bb_xyz_p1']=xyz[:,:3]
        else:
            feature_dict['sequence_p1']=seq[:self.binderlen]
            feature_dict['sequence_p2']=seq[self.binderlen:]
            feature_dict['bb_xyz_p1']=xyz[:self.binderlen,:3]
            feature_dict['bb_xyz_p2']=xyz[self.binderlen:,:3]

        return self.recover_sc(feature_dict)

    
    def __call__(self, xyz, seq=None, seq_mask=None, ligand_xyz=None, ligand_aatypes=None):

        xyz=xyz.clone().to(self.device)
        
        if seq!=None:
            seq=seq.clone().to(self.device)
            seq_mask=seq_mask.clone().to(self.device)

        if xyz.shape[0]<=self.binderlen:
            self.binderlen=-1
        
        if self.recover_sc==None:
            self.init_recover_sc()

        seq=self.get_aa_probs(xyz, seq, seq_mask, ligand_xyz, ligand_aatypes)

        d=self.bb2martini(xyz, seq)

        return d


class BB2MartiniModule(nn.Module):
    
    def __init__(self, chains=['_p1','_p2'], molecule='protein', device='cpu'):

        nn.Module.__init__(self)
        
        self.chains=chains
        
        # get ideal pseudoatom coords and types for aminoacids
        if molecule=='protein':
            ps, ts=martinize(num2aa, 
                     [list(x.keys()) for x in ideal_coords],
                     [list(x.values()) for x in ideal_coords] )
        if molecule=='na':
            ps, ts=martinize(num2na, 
                     [list(x.keys()) for x in ideal_coords_na],
                     [list(x.values()) for x in ideal_coords_na] )
        assert len(ps)==len(ts)

        self.num_aa=len(ts) # A
        self.num_pseudo=max([len(x) for x in ts]) # P

        # encode types
        num2ps=list(set([x for y in ts for x in y ]))
        ps2num={x:i for i, x in enumerate(num2ps)}

        rs=[[pseudoatom_radii[y] for y in x] for x in ts]
        ws=[[pseudoatom_weights[y] for y in x] for x in ts]
        ts=[[ps2num[y] for y in x] for x in ts] 

        self.num_types=len(num2ps) # T
        
        # pad coords and types
        for i in range(len(ps)):
            assert len(ps[i])==len(ts[i])
            to_add=self.num_pseudo-len(ps[i])
            for j in range(to_add):
                ps[i]=[ps[i][0]]+ps[i]
                ts[i]=[ts[i][0]]+ts[i]
                rs[i]=[rs[i][0]]+rs[i]
                ws[i]=[ws[i][0]]+ws[i]
            for j in range(to_add):
                ws[i][j]/=(to_add+1)

        ps=np.array(ps)
        ts=np.array(ts)
        rs=np.array(rs)
        ws=np.array(ws)
        
        self.map_coords=torch.Tensor(ps)[None,:,:,:].to(device) # (1,A,P,3)
        self.map_types=F.one_hot(torch.Tensor(ts).to(torch.int64),
                                 num_classes=self.num_types).float()[None,:,:,:].to(device) # (1,A,P,T)
        self.map_weights=torch.Tensor(ws)[None,:,:,None].to(device) # (1,A,P,1)
        self.map_radii=torch.Tensor(rs)[None,:,:,None].to(device) # (1,A,P,1)

        self.map_coords=nn.Parameter(self.map_coords, requires_grad=False)
        self.map_types=nn.Parameter(self.map_types, requires_grad=False)
        self.map_weights=nn.Parameter(self.map_weights, requires_grad=False)
        self.map_radii=nn.Parameter(self.map_radii, requires_grad=False)

    def forward(self, data):

        for ch in self.chains: 
            sequence=data[f'sequence{ch}'][:,:,None,None]
            bb=data[f'bb_xyz{ch}']
                       
            xyz = (sequence*self.map_coords*self.map_weights).sum(1)/(sequence*self.map_weights).sum(1) # (L,P, 3)      
            types = (sequence*self.map_types*self.map_weights).sum(1)/(sequence*self.map_weights).sum(1) # (L, P, T)
            radii = (sequence*self.map_radii*self.map_weights).sum(1)/(sequence*self.map_weights).sum(1) # (L, P, 1)
            weights = (sequence*self.map_weights).sum(1) # (L, P, 1)
        
            Rs, Ts = rigid_from_3_points(bb[None,:,0,:],bb[None,:,1,:],bb[None,:,2,:]) # transforms
            xyz = torch.einsum('lpij,lpj->lpi', Rs, xyz.transpose(0,1)) + Ts
        
            xyz=xyz.transpose(0,1).reshape((-1,3))
            types=types.reshape((-1,self.num_types))
            radii=radii.reshape(-1)
            weights=weights.reshape(-1)

            data[f'atom_xyz{ch}']=xyz
            data[f'atom_types{ch}']=types
            data[f'atom_rad{ch}']=radii
            data[f'atom_weights{ch}']=weights

        return data 


class ReshapeBB: # то сut backbone atoms from the protein

    def __init__(self, chains=['_p1','_p2']):
        self.chains=chains
    
    @torch.no_grad()
    def __call__(self,data):

        for ch in self.chains:                      
            seq=torch.stack([data[f'atom_chain{ch}'],data[f'atom_resid{ch}'], data[f'atom_type{ch}']], dim=1)
            seq, idx=seq.unique(dim=0, sorted=True, return_inverse=True)
            idx=(idx[:,None]==torch.arange(len(seq),device=idx.device)[None,:]).to(int).argmax(dim=0)         
            for key in data.keys:
                if ch not in key:
                    continue
                data[key]=data[key][idx]
            
            seq1=seq[:,:2].unique(dim=0)    
            if len(seq1)*3!=len(data[f'atom_resid{ch}']):
                mask=torch.ones_like(data[f'atom_chain{ch}']).to(bool)
                for aa in seq:
                    idx=(data[f'atom_resid{ch}']==aa[1])*(data[f'atom_chain{ch}']==aa[0])
                    if sum(idx)<3:
                        mask=~idx*mask
                for key in data.keys:
                    if ch not in key:
                        continue
                    data[key]=data[key][mask]   
                
                seq=seq[mask,:2].unique(dim=0)
                assert len(seq)*3==len(data[f'atom_resid{ch}'])
    
            bb_xyz=torch.stack((data[f'atom_xyz{ch}'][0::3],
                                data[f'atom_xyz{ch}'][1::3],
                                data[f'atom_xyz{ch}'][2::3]), 
                               dim=1)
            for key in list(data.keys):
                if ch not in key:
                    continue
                elif 'atom' in key:
                    del data[key]
                else:
                    data[key]=data[key][::3].detach()
            data[f'bb_xyz{ch}']=bb_xyz.detach() # coordinates of 3 backbone atoms

        return data    
    