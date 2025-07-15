import os
import sys
from subprocess import run, CalledProcessError

class SubmitTrajectories:

    def __init__(self, partition):
        self.partition = partition

    def read_sampling(self):
        prop = open("sampling.inp", 'r+')
        for line in prop:
            if "model = " in line:
                model = str(line.split()[2])
            elif "from = molden" in line:
                model = "molecule"
        return model 

    def read_spp(self):
        prop = open("spp.inp", 'r+')
        for line in prop:
            if "use_db =" in line:
                use_db = str(line.split()[2])
        return use_db 

    def read_prop(self):
        prop = open("prop.inp", 'r+')
        for line in prop:
            if "method =" in line:
                method = str(line.split()[2])
        return method 
    
    def get_best_partition(self):
        
        excluded_partitions = {'eiguren', 'bcam-exclusive', 'biophys', 'garciaetxarri', 'chen'}
        preferred_states = ['idle', 'mix']
        partitions = []
        
        output = subprocess.check_output(['sinfo', '-o', '%P %D %t %N'],
        text=True)

        for line in output.splitlines()[1:]:  # Skip header
            parts = line.split()
            if len(parts) < 4:
                continue

            raw_partition, nodes_str, state, _ = parts
            partition = raw_partition.rstrip('*')
            nodes = int(nodes_str)

            if partition in excluded_partitions:
                continue

            if state.lower() in preferred_states:
                # Lower index in preferred_states means higher priority (idle=0, mix=1)
                partitions.append((preferred_states.index(state.lower()),-nodes, partition))
        if partitions:
            partitions.sort()
            return partitions[0][2]  # Best match

        return "general"  # Fallback


    def trajectories(self):
        model = self.read_sampling()
        use_db = self.read_spp()
        method = self.read_prop()
        for traj in os.listdir("prop"):
            subfolder = os.path.join("prop",traj)
            partition = self.get_best_partition()
            print(f"Using partition: {partition}")
            try:
                if method == "Born_Oppenheimer":
                    run(['sbatch',
                    f'--partition={partition}',
                    'dipc_atlas_run.sh'], cwd=subfolder, check=True, shell=True)
                elif method == "LandauZener":
                    if use_db == "yes":
                        run(['sbatch db_lz_run.sh'], cwd=subfolder, check=True, shell=True)
                    else:
                        run(['sbatch lz_run_om.sh'], cwd=subfolder, check=True, shell=True)
                        #run(['sbatch lz_run.sh'], cwd=subfolder, check=True, shell=True)
                elif method == "Surface_Hopping":
                    if model == "Tully_1":
                        run(['sbatch model_fssh_run.sh'], cwd=subfolder, check=True, shell=True)
                    elif model == "LVC":
                        run(['sbatch lvc_fssh_run.sh'], cwd=subfolder, check=True, shell=True)
                    #else:
                    #    #run(['sbatch lz_nacs.sh'], cwd=subfolder, check=True, shell=True)
                    #    #run(['sbatch saoovqe.sh'], cwd=subfolder, check=True, shell=True)
                    #    if self.partition == "cpu_lorentz_omolcas":
                    #        run(['sbatch cpu_lorentz_open_molcas.sh'], cwd=subfolder, check=True, shell=True)
                    #    elif self.partition == "cpu_long_vqe":
                    #        run(['sbatch cpu_long_saoovqe.sh'], cwd=subfolder, check=True, shell=True)
                    #    elif self.partition == "cpu_lorentz_vqe":
                    #        run(['sbatch cpu_lorentz_saoovqe.sh'], cwd=subfolder, check=True, shell=True)
                    #    #run(['sbatch test_open_molcas.sh'], cwd=subfolder, check=True, shell=True)
                    #    #run(['sbatch open_molcas.sh'], cwd=subfolder, check=True, shell=True)
                    #    #run(['sbatch bagel.sh'], cwd=subfolder, check=True, shell=True)
                    #    #run(['sbatch lvc_fssh_run.sh'], cwd=subfolder, check=True, shell=True)
                    #    #run(['sbatch fssh_run.sh'], cwd=subfolder, check=True, shell=True)
            except KeyboardInterrupt or CalledProcessError:
                break
            print("Submitting", subfolder)            


if __name__=='__main__':
    partition = sys.argv[1]
    all_traj = SubmitTrajectories(partition)
    all_traj.trajectories()
