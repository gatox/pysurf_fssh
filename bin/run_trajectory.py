from colt import from_commandline
from pysurf.fssh import State, VelocityVerletPropagator
from pathlib import Path


import os
path_test = "/Users/salazar/Desktop/new_test_nofvqe/h2_noisless/prop/traj_test"
os.chdir(path)


@from_commandline("""
inputfile = prop.inp :: file
""")
def command_run_trajectory(inputfile="prop.inp"):
    db_file = "results.db"
    restart = Path(db_file).exists() #and Path("gen_results.out").exists() 
    if restart:
        print("Restarting trajectory from last md_step ...")
        elec_state = State.from_db_frame(db_file, config_file=inputfile)
    else:
        elec_state = State.from_questions(config=inputfile)
    DY = VelocityVerletPropagator(elec_state, restart=restart)
    try:
        result_2 = DY.run()
    except SystemExit as err:
        print("An error:", err)

if __name__=="__main__":
    command_run_trajectory()
