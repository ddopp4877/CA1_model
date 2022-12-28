# Model of CA1 section of Hippocampus
### Goal is to show SPWR in CA1
to make an environment where everything will run:
python3 -m venv CA1_model_env
source CA1_model_env/bin/activate
pip install -r requirements.txt

## Building network
```
python3 build_network_V2.py
if there is an issue with output then just
rm -rf output
mkdir output
to build
alternatively can do all this with:
bash eServerBuild.sh
```

### Running network
```
python3 run_network.py


currently three AACs are stimulated
background noise is currently not used.



using build_env_bionet will create a circuit json which the from_config in run_network.py cannot find the modfiles for some reason.
but, the build_env_bionet() will create a proper simulation_config.json, so run build_network as is,
then add mechanisms_dir": "$COMPONENTS_DIR\\mechanisms\\modfiles",
to the circuit_config.json file
