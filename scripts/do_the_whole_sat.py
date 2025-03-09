import os
import subprocess

benchmark = ["b19", "b20", "b21", "b22"]

def execute_program(input_files, output_folder, folder_name, exe_folder):
    # 构建执行命令
    command = os.path.join(exe_folder, "TestATPG ")
    command += input_files + " > " + os.path.join(exe_folder, folder_name + ".report")

    # 执行命令
    subprocess.run(command, shell=True)

def main():
    root_folder = "/home/dell/Desktop/test_circuit2/"
    exe_folder = "../cmake-build-debug/test/"

    # 遍历文件夹
    for folder_name in benchmark:
        folder_path = os.path.join(root_folder, folder_name)
        if os.path.isdir(folder_path) and folder_name.startswith("b"):
            # 构建输入文件路径
            input_files = str(os.path.join(folder_path, folder_name + "_opt_C_primitives.vy")) + " " + \
                          str(os.path.join(folder_path, folder_name + "_opt_C_cells.vy")) + " " + \
                          str(os.path.join(folder_path, folder_name + "_opt_C.spf")) +  " " + \
                          folder_path
            output_folder = os.path.join(root_folder, folder_name)
            # 确保输出文件夹存在
            os.makedirs(output_folder, exist_ok=True)
            # 执行可执行程序
            execute_program(input_files, output_folder, folder_name, exe_folder)

if __name__ == "__main__":
    main()
