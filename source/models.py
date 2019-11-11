import torch
from source.parameters import dbsearch
from source.parameters import types

# Define types of variables from parameters
dtype = types['dtype']
torch_device = types['torch_device']
float_type = types['float_type']
# Define log-softmax function in advance
log_softmax = torch.nn.LogSoftmax(dim=0)

class RBMModel(torch.nn.Module):
    """
    Make an RBM model
    """
    def __init__(self, dbs):
    
        super(RBMModel, self).__init__()
        self.dbs = dbs
        self.rbm_dim = self.dbs.max_bin
        self.rbm_weights = torch.eye(self.rbm_dim,   dtype=float_type, device=torch_device, requires_grad=True)
        self.rbm_vb      = torch.zeros(self.rbm_dim, dtype=float_type, device=torch_device, requires_grad=True)
        self.rbm_hb      = torch.zeros(self.rbm_dim, dtype=float_type, device=torch_device, requires_grad=True)
        self.diagonal_mask = torch.eye(self.rbm_dim, dtype=torch.uint8, device=torch_device, requires_grad=False)

    def state_dict(self):
        """
        Provides ability to save only state_dict instead of whole the model

        Notice: analog of built-in method torch.save(model.state_dict())
        """
        return {"weight" : self.rbm_weights, "vb" : self.rbm_vb, "hb" : self.rbm_hb, "dim" : self.rbm_dim}
    
    def weight_update(self, loss, optimizer, lr=0.001):
        """Updating weights with customized parameters and preserving masked diagonal matrix"""
        params = [self.rbm_weights, self.rbm_vb, self.rbm_hb]
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        self.rbm_weights.data[self.diagonal_mask] = 1
        
    def scoring(self, experimental_spectra, theoretical_spectra, candidate_spectra):
        """
        Applies vectorized scoring with positive and negative phases
        
        Notice: using candidate_spectra instead of performing MCMC sampling
        """
        visible_data = torch.from_numpy(experimental_spectra).type(dtype).to(torch_device)
        candidate_data = torch.from_numpy(candidate_spectra).type(dtype).to(torch_device)

        transformed = torch.mm(visible_data, self.rbm_weights)
        visible_bias = torch.mv(visible_data, self.rbm_vb)
        candidate_bias = torch.mv(candidate_data, self.rbm_hb)
        
        scores = []

        for i, theo_data in enumerate(theoretical_spectra):
            theo_pept_num = len(theo_data)
            if theo_pept_num == 0:
                continue
            positive_phase = torch.zeros((1), dtype=float_type).to(torch_device)
            
            positive_phase += torch.dot(transformed[i], candidate_data[i]) \
                            + visible_bias[i] \
                            + candidate_bias[i]

            negative_phase = torch.mv(theo_data, transformed[i]) \
                            + torch.mv(theo_data, self.rbm_hb) \
                            + visible_bias[i]

            scores.append(torch.cat((positive_phase, negative_phase)))
        
        return scores
