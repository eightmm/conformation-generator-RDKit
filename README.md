# Conformation Generator

A Python tool for generating molecular conformations using RDKit's ETKDGv3 algorithm and optimizing them with UFF and MMFF force fields.

This project provides a simple and flexible way to generate multiple conformers of molecules, align them, optimize with different force fields, and output the results in SDF format with RMSD values for each conformer.

## Table of Contents
1. [Installation](#installation)
2. [Usage](#usage)
3. [Command Line Arguments](#command-line-arguments)
4. [Example](#example)
5. [Output](#output)
6. [License](#license)

---

## Installation

1. First, clone the repository:
    \`\`\`bash
    git clone https://github.com/your-username/conformation-generator.git
    cd conformation-generator
    \`\`\`

2. Install the required dependencies (RDKit):
    RDKit can be installed via conda:
    \`\`\`bash
    conda create -c conda-forge -n rdkit-env rdkit
    conda activate rdkit-env
    \`\`\`

3. If you're using another environment, ensure that RDKit and Python dependencies are properly installed.

---

## Usage

You can generate molecular conformations by running the \`generate_conformer.py\` script. 

### Basic Usage:
\`\`\`bash
python generate_conformer.py --sdf <input_sdf> --out <output_sdf> --num <number_of_conformers> --max_iter <max_iterations> --n_cpu <number_of_cpu> --verbose
\`\`\`

### Example:
\`\`\`bash
python generate_conformer.py --sdf example/molecule.sdf --out example/molecule_conf.sdf --num 100 --max_iter 200 --n_cpu 4 --verbose
\`\`\`

### Detailed Steps:
1. **Input**: The script accepts a molecule in SDF format (\`.sdf\`) as input.
2. **Conformer Generation**: It generates multiple conformers using the ETKDGv3 algorithm.
3. **Optimization**: Two force fields (UFF and MMFF) are applied to optimize the generated conformers.
4. **Alignment and RMSD Calculation**: The conformers are aligned, and the RMSD values for each conformer are computed.
5. **Output**: The results are saved in an SDF file with each conformer and its associated RMSD value.

---

## Command Line Arguments

| Argument        | Description                                                                 | Default         | Required |
|-----------------|-----------------------------------------------------------------------------|-----------------|----------|
| \`--sdf\`         | Input SDF file with the molecule                                             | N/A             | Yes      |
| \`--out\`         | Output SDF file to save the generated conformers                             | \`<input>_conf.sdf\` | No    |
| \`--num\`         | Number of conformations to generate                                          | 100             | No       |
| \`--max_iter\`    | Maximum iterations for force field optimization                              | 200             | No       |
| \`--n_cpu\`       | Number of CPU cores to use for the calculations                              | 0 (all cores)   | No       |
| \`--verbose\`     | Enable verbose output to display detailed information during execution       | False           | No       |

---

## Example

This repository comes with an example directory that includes a sample molecule (\`molecule.sdf\`) and a sample run script (\`run.sh\`). 

To try it out, simply run:

\`\`\`bash
cd example
bash run.sh
\`\`\`

Here’s the breakdown of what happens:
- The script reads the input molecule from \`molecule.sdf\`.
- It generates 100 conformers of the molecule.
- UFF and MMFF optimizations are applied.
- Conformers are aligned and written to \`molecule_conf.sdf\`, with the RMSD values stored in each conformer entry.

---

## Output

The output will be an SDF file containing:
- All the generated conformers.
- RMSD values for each conformer, calculated based on the alignment with the reference conformer.

Each conformer in the output SDF file will have a property named \`RMSD\` that indicates how far it is (in terms of RMSD) from the first conformer. This allows you to easily assess the diversity of the generated conformations.

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## Contributing

Feel free to open issues or submit pull requests if you encounter any problems or have suggestions for improvements.

