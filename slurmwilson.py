import paramiko
import os

def get_slurm_token():
    user = os.getenv("USER")
    client = paramiko.SSHClient()
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    client.connect("wilson")
    command = "scontrol token lifespan=259200" # is this unreasonable?
    stdin, stdout, stderr = client.exec_command(command)
    output = stdout.read().decode("utf-8")
    token = None
    for line in output.split("\n"):
        if line.startswith("SLURM_JWT"):
            token = line.split("=")[1].strip()
        else:
            pass
    client.close()
    return user, token

user, token = get_slurm_token()