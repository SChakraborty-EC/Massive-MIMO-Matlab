function [BER, node_count]=plot4_Compile(plot_data, number)

consolidated_data  = consolidate(plot_data,number);

% save the consolidate data into a matrix for export

[BER, node_count] = data_save(consolidated_data);

end

function [BER, node_count] = data_save(result)
BER = zeros(length(result.free_nodes), 2);
BER(:,1) = result.free_nodes;
BER(:,2) = result.BER;

node_count = zeros(length(result.free_nodes), 2);
node_count(:,1) = result.free_nodes;
node_count(:,2) = result.count;
end


% =========================================================================
% Title       : Consolidate Simulation results
% File        : consolidate_BER.m 
% =========================================================================

function Out = consolidate(filename,number)
  disp(sprintf('consolidate %s ...',filename));
  k = 0;
  for i=1:number
    fullname = ['Results/',filename,'_',num2str(i-1),'.mat'];
    % -- check whether file exists
    if exist(fullname)>0
      tmp = load(fullname);
      if k==0
        ptype = tmp;
        BER = ptype.Results.BER;
        count = ptype.Results.node_count;
      else
        count = count + tmp.Results.node_count;
      end
      k=k+1;
    end
  
  end
  % -- average and compute output
  ptype.Results.BER = BER/k;
  ptype.Results.count = count/k;
  Out = ptype.Results;
end
