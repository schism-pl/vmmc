#!/bin/bash
  
#Setting number of toml files
#Note that "start_radius" is the next radius increment than the one specified in "h00.toml" below
max_num=36
start_radius=0.02
end_radius=0.37
increment=0.01

proj_dir="{dir_name_here}"
echo "removing previous directory $proj_dir"
rm -r $proj_dir
mkdir -p $proj_dir

#Creating Binding Energy Folders
for i in $(seq -w 1 20) #Sets the range of binding energies to be explored (1-20 in this case)
do
        dir_path="$proj_dir/be$i"
        echo "Creating input directory $dir_path"
        mkdir -p $dir_path
        cat > "$dir_path/h00.toml" << EOF 
seed = 168235197561895567
initial_particles = 500
box_width = 30.0
box_height = 30.0

[protocol]
chemical_potential_eq = "0.0"
interaction_energy_eq = "$i.0"
num_megasteps = 1000

[[shapes]]

[[shapes.patches]]
radius = 0.01
theta = 0.0
chemtype = 0

[[shapes.patches]]
radius = 0.01
theta = 60.0
chemtype = 0

[[shapes.patches]]
radius = 0.01
theta = 120.0
chemtype = 0

[[shapes.patches]]
radius = 0.01
theta = 180.0
chemtype = 0

[[shapes.patches]]
radius = 0.01
theta = 240.0
chemtype = 0

[[shapes.patches]]
radius = 0.01
theta = 300.0
chemtype = 0
EOF
done

echo "Creating toml files"

for i in $(seq -w 1 20) 
do
        dir_path="$proj_dir/be$i"
        current_radius=$start_radius 
        for j in $(seq -w 0 $max_num) 
        do
                if [ "$j" != "00" ]; then
                        cp "$dir_path/h00.toml" "$dir_path/h$j.toml"
                        formatted_radius=$(printf "%.2f" $current_radius)
                        sed -i "s|radius = 0.01|radius = $formatted_radius|" "$dir_path/h$j.toml"
                        current_radius=$(echo "$current_radius + $increment" | bc) #check if you have bc installed
                fi 
                echo "running $dir_path/h$j.toml"
                /{path_to_vmmc_executable}/vmmc/target/release/vmmc-run --input="$dir_path/h$j.toml" --output-dir="$dir_path/h$j" &> "$dir_path/h$j.log" &
        done
        echo "Finished creating $dir_path"
done
echo "DONE running vmmc in $dir_path"

wait


