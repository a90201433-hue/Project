for cfg in config_set/*.toml; do
    echo "Making pictures for $cfg"

    name=$(basename "$cfg" .toml)
    python3 PlotSlice.py "$name"
    for field in speed P rho u v; do
        python3 PlotMap.py "$name" "$field"
    done 
done
