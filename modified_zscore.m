function m_z = modified_zscore(x)
m_z = 0.6745*(x-median(x))/mad(x,1);
end