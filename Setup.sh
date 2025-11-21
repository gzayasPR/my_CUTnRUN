
proj_dir="$(pwd)"
my_bin="$(pwd)/bin"
my_data="$(pwd)/data"
my_results="$(pwd)/results"
my_docs="$(pwd)/docs"
my_bash="$(pwd)/bash_out"
my_softwares="$(pwd)/softwares"
# Create directories
mkdir -p "$my_bin" "$my_data" "$my_results" "$my_docs" "$my_bash" "$my_softwares"

env_file="$my_bin/project_env.sh"


# Write project environment file
env_file="${my_bin}/project_env.sh"
cat <<EOL > "$env_file"
#!/bin/bash
export proj_dir="${proj_dir}"
export my_bin="${my_bin}"
export my_data="${my_data}"
export my_results="${my_results}"
export my_docs="${my_docs}"
export my_bash="${my_bash}"
export my_softwares="${my_softwares}"
EOL

# Make the file executable
chmod +x "$env_file"


cd $my_bin
