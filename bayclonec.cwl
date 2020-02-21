cwlVersion: v1.0
class: CommandLineTool
label: BayCloneC
baseCommand: ["Rscript", "./run_bc.sh"]
requirements:
  - class: DockerRequirement
    dockerPull: subhajit0606/bayclonec:v1

inputs:
  input_vcf:
    type: File
    inputBinding:
      position: 1

  battenberg_file:
    type: File
    inputBinding:
      position: 2

  purity_file:
    type: File
    inputBinding:
      position: 3

  prefix:
    type: string
    inputBinding:
      position: 4

outputs:
  population:
    type: File
    outputBinding:
      glob: 1B.txt

  proportion:
    type: File
    outputBinding:
      glob: 1C.txt

  cluster_assignment:
    type: File
    outputBinding:
      glob: 2A.txt

  cocluster_assignment:
    type: File
    outputBinding:
      glob: 2B.txt
