#!/bin/bash

# Collect source data used in Donovan et al., 2024 (https://www.nature.com/articles/s41467-024-49781-1)


cd "/Users/jessica/Dropbox/EspinosaGroup/ANALYSIS/2025/dev/rawdata";

synapse login "shawjr" "Bloombug2015!";

synapse get -q "SELECT * FROM syn31481952.8";

echo "done"