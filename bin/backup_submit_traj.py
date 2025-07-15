import os
import sys
import subprocess
from subprocess import run, CalledProcessError

class SubmitTrajectories:

    def __init__(self, partition):
        self.partition = partition

    def read_sampling(self):
        with open("sampling.inp", 'r') as prop:
            for line in prop:
                if "model = " in line:
                    model = str(line.split()[2])
                elif "from = molden" in line:
                    model = "molecule"
        return model

    def read_spp(self):
        with open("spp.inp", 'r') as prop:
            for line in prop:
                if "use_db =" in line:
                    use_db = str(line.split()[2])
        return use_db

    def read_prop(self):
        with open("prop.inp", 'r') as prop:
            for line in prop:
                if "method =" in line:
                    method = str(line.split()[2])
        return method

    def get_best_partition(self):
        excluded = {'eiguren', 'bcam-exclusive', 'biophys', 'garciaetxarri', 'chen'}
        preferred_partitions = ['cpu-long', 'cpu_lorentz', 'general']
        preferred_states = ['idle', 'mix']
        partitions = []

        try:
            output = subprocess.check_output(['sinfo', '-o', '%P %D %T %N %l'], text=True)
        except subprocess.CalledProcessError:
            print("Error running sinfo")
            return "general", "7-00:00:00"

        for line in output.splitlines()[1:]:
            parts = line.split()
            if len(parts) < 5:
                continue

            raw_partition, nodes_str, state, _, time_limit = parts
            partition = raw_partition.rstrip('*')
            if partition in excluded:
                continue
            if state.lower() in preferred_states:
                partitions.append((
                    preferred_states.index(state.lower()),
                    -int(nodes_str),
                    preferred_partitions.index(partition) if partition in preferred_partitions else 99,
                    partition,
                    time_limit
                ))

        if partitions:
            partitions.sort()
            _, _, _, best_partition, time_limit = partitions[0]
            # Sanitize and validate time limit
            if time_limit in ["infinite", "n/a", "UNKNOWN", "None", None]:
                time_limit = "7-00:00:00"  # Default max if unspecified or invalid
            return best_partition, time_limit

        return "general", "7-00:00:00"

    def trajectories(self):
        model = self.read_sampling()
        use_db = self.read_spp()
        method = self.read_prop()

        for traj in os.listdir("prop"):
            subfolder = os.path.join("prop", traj)
            #partition, time_limit = self.get_best_partition()
            #partition, time_limit = "general", "7-00:00:00" 
            partition, time_limit = "regular", "1-00:00:00" 
            print(f"Using partition: {partition} with time limit: {time_limit}")

            try:
                if method == "Born_Oppenheimer":
                    run(['sbatch',
                        #f'--partition={partition}',
                        #f'--time={time_limit}',
                        'pysurf_pynof_run.sh'
                    ], cwd=subfolder, check=True,shell=True)

                elif method == "LandauZener":
                    if use_db == "yes":
                        run(['sbatch', 'db_lz_run.sh'], cwd=subfolder, check=True)
                    else:
                        run(['sbatch', 'lz_run_om.sh'], cwd=subfolder, check=True)

                elif method == "Surface_Hopping":
                    if model == "Tully_1":
                        run(['sbatch', 'model_fssh_run.sh'], cwd=subfolder, check=True)
                    elif model == "LVC":
                        run(['sbatch', 'lvc_fssh_run.sh'], cwd=subfolder, check=True)
                    else:
                        #run(['sbatch lz_nacs.sh'], cwd=subfolder, check=True, shell=True)
                        #run(['sbatch saoovqe.sh'], cwd=subfolder, check=True, shell=True)
                        if self.partition == "cpu_lorentz_omolcas":
                            run(['sbatch cpu_lorentz_open_molcas.sh'], cwd=subfolder, check=True, shell=True)
                        elif self.partition == "cpu_long_vqe":
                            run(['sbatch cpu_long_saoovqe.sh'], cwd=subfolder, check=True, shell=True)
                        elif self.partition == "cpu_lorentz_vqe":
                            run(['sbatch cpu_lorentz_saoovqe.sh'], cwd=subfolder, check=True, shell=True)
                        #run(['sbatch test_open_molcas.sh'], cwd=subfolder, check=True, shell=True)
                        #run(['sbatch open_molcas.sh'], cwd=subfolder, check=True, shell=True)
                        #run(['sbatch bagel.sh'], cwd=subfolder, check=True, shell=True)
                        #run(['sbatch lvc_fssh_run.sh'], cwd=subfolder, check=True, shell=True)
                        #run(['sbatch fssh_run.sh'], cwd=subfolder, check=True, shell=True)
            except (KeyboardInterrupt, CalledProcessError) as e:
                print(f"Submission interrupted or failed in {subfolder}: {e}")
                break

            print("Submitting", subfolder)


if __name__ == '__main__':
    partition = sys.argv[1] if len(sys.argv) > 1 else None
    all_traj = SubmitTrajectories(partition)
    all_traj.trajectories()
