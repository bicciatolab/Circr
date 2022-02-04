### Run on server with Ubuntu 18.04

### Circr folder downloaded and saved in home folder

# running Circr.py using the full sequence coordinates
python3 ~/Circr/Circr.py -i circular_mixed.bed -s mouse -v mm9 -o test_circr_server.txt

# running Circr.py using the final features coordinates
python3 ~/Circr/Circr.py -i circular_mixed.bed -c -s mouse -v mm9 -o test_circr_server_coord.txt

