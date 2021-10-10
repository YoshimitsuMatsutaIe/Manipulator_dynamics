# import time
# for i in range(50):
#     print("\r[{}]".format("#"*i), end="")
#     time.sleep(0.1)


import time
epoch = 50
for i in range(epoch):
    bar = "="*i + (">" if i < epoch-1 else "=") + " "*(epoch-i-1)
    print("\r[{}] {}/{}".format(bar, i+1, epoch), end="")
    #time.sleep(0.1)