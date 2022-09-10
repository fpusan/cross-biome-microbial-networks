00-algorithm/core-linkage.py -q Thermal -r 5 -m 5 -s 10 -c 0.0001 -b 100 -p subtype -k 0 -g 10 -y 6 -z 5 -o thermal
	#Note the -r 5 for extreme environments
00-algorithm/core-linkage.py -q Saline -w Marine -e Sediment -r 5 -m 5 -s 10 -c 0.0001 -b 100 -p subtype -k 0 -g 10 -y 6 -z 5 -o marine-sediment
00-algorithm/core-linkage.py -q Saline -w Marine -e Water -r 5 -m 5 -s 10 -c 0.0001 -b 100 -p subtype -k 0 -g 10 -y 6 -z 5 -o marine-water
00-algorithm/core-linkage.py -q Non-saline -w Soils -r 5 -m 5 -s 10 -c 0.0001 -b 100 -p subtype -k 0 -g 10 -y 6 -z 5 -o soils
#00-algorithm/core-linkage.py -q Host-associated -w Gut -r 10 -m 10 -s 10 -c 0.0001 -b 100 -p subtype -k 0 -g 10 -y 6 -z 5 -o host-gut
#00-algorithm/core-linkage.py -q Host-associated -w Oral -r 10 -m 10 -s 10 -c 0.0001 -b 100 -p subtype -k 0 -g 10 -y 6 -z 5 -o host-oral
00-algorithm/core-linkage.py -q Artificial -w "Water treatment" -r 5 -m 5 -s 10 -c 0.0001 -b 100 -p subtype -k 0 -g 10 -y 6 -z 5 -o water-treatment
00-algorithm/core-linkage.py -q Saline -w Hypersaline -r 5 -m 5 -s 10 -c 0.0001 -b 100 -p subtype -k 0 -g 10 -y 6 -z 5 -o saline-hypersaline
	#Lower richness cutoff for extreme environment
00-algorithm/core-linkage.py -q Hypothermal -w Non-marine -r 5 -m 5 -s 10 -c 0.0001 -b 100 -p subtype -k 0 -g 10 -y 6 -z 5 -o hypothermal-nonmarine
        #Lower richness cutoff for extreme environment
00-algorithm/core-linkage.py -q Other -w Oil -r 5 -m 5 -s 10 -c 0.0001 -b 100 -p subtype -k 0 -g 10 -y 6 -z 5 -o oil
        #Lower richness cutoff for extreme environment
00-algorithm/core-linkage.py -q Non-saline -w Freshwaters -r 5 -m 5 -s 10 -c 0.0001 -b 100 -p subtype -k 0 -g 10 -y 6 -z 5 -o freshwater
00-algorithm/core-linkage.py -q Host-associated -r 5 -m 5 -s 10 -c 0.0001 -b 100 -p subtype -k 0 -g 10 -y 6 -z 5 -o host-associated
	#Two null models fail to fit
