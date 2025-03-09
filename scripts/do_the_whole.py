import os
import subprocess

def execute_program(input_files, output_folder, folder_name):
    # 构建执行命令
    command = "../cmake-build-debug/test/TestATPG "
    command += input_files + " > " + os.path.join(output_folder, folder_name + ".report")

    # 执行命令
    subprocess.run(command, shell=True)

def main():
    root_folder = "../example/data/new_bench/"


    # 遍历文件夹
    for folder_name in os.listdir(root_folder):
        folder_path = os.path.join(root_folder, folder_name)
        if os.path.isdir(folder_path) and (folder_name.startswith("s") or folder_name.startswith("c")):
            # 构建输入文件路径
            input_files = str(os.path.join(folder_path, folder_name + "_primitives.vy")) + " " + \
                          str(os.path.join(folder_path, folder_name + "_cells.vy")) + " " + \
                              str(os.path.join(folder_path, folder_name + ".spf")) +  " " + \
                                  str(os.path.join(folder_path, "out"))
            output_folder = os.path.join(root_folder, folder_name, "out")
            # 确保输出文件夹存在
            os.makedirs(output_folder, exist_ok=True)
            # 执行可执行程序
            execute_program(input_files, output_folder, folder_name)

if __name__ == "__main__":
    main()
