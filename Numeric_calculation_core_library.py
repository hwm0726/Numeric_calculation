#------------------------------------------------------------------------------
#算法1：列选主元消去法解线性方程组

def column_selected(a,n,x):  #参数依次为：增广矩阵a,阶数n,一维数组x（存放结果）
    #消元n次
    for i in range(n):
        #列选主元
        k = i
        for j in range(i,n):
            if abs(a[j][i]) > abs(a[k][i]):
                k = j
        #k行与i行互换
        tmp = 0
        if k != i:    
            for j in range(i,n + 1):
                tmp = a[i][j]
                a[i][j] = a[k][j]
                a[k][j] = tmp
        #i行归一化
        tmp = a[i][i]
        for j in range(i,n + 1):
            a[i][j] = a[i][j] / tmp
        #(i+1)行至第n行消元
        for j in range(i + 1,n):
            tmp = a[j][i]
            for k in range(i,n + 1):
                a[j][k] = a[j][k] - tmp * a[i][k]
    #逆序回代解出结果，存至数组x
    for i in range(n - 1,-1,-1):
        x[i] = a[i][n] / a[i][i]
        if i > 0:
            for j in range(n - 1,i - 1,-1):
                a[i - 1][n] = a[i - 1][n] - a[i - 1][j] * x[j]
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
#算法2：Doolitter（杜立特尔）分解

def doolitter(a,n,l,u):#参数依次为：方阵a,阶数n,二维矩阵l和u（存放结果）
    #计算u的第一行,l的第一列
    for j in range(n):
        u[0][j] = a[0][j]
    for i in range(1,n):
        l[i][0] = a[i][0] / u[0][0]
    for r in range(1,n):
        #计算u第r行
        for j in range(r,n):
            u[r][j] = a[r][j]
            for k in range(r):
                u[r][j] = u[r][j] - l[r][k] * u[k][j]
        #计算l第r列
        for i in range(r + 1,n):
            l[i][r] = a[i][r]
            for k in range(r):
                l[i][r] = l[i][r] - l[i][k] * u[k][r]
            l[i][r] = l[i][r] / u[r][r]
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
#算法3：拉格朗日插值法

def langrange(x,y,n,xi):#参数依次为：矩阵x,y,插值点个数n,待求点的x值为xi
    #res为结果
    res = 0.0
    #代入拉格朗日公式
    for k in range(n):
        #求第k个插值函数
        l = 1.0
        for j in range(n):
            if j != k:
                 l = l * (xi - x[j]) / (x[k] - x[j])
        #累加
        res = res + l * y[k]
    #返回计算结果
    return res
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
#算法4：牛顿插值法

def newton(x,y,n,xi):#参数依次为：矩阵x,y,插值点个数n,待求点的x值为xi
    #f为一维差商数组
    f = []
    for i in range(n):
        f.append(y[i])
    #计算1~n阶差商
    for i in range(1,n):
        #逆序计算差商
        for j in range(n - 1,i - 1,-1):
            f[j] = (f[j - 1] - f[j]) / (x[j - i] - x[j])
    #结果累加器
    res = 0.0
    mul = 1.0 #秦九韶迭代器
    for i in range(n):
        if i >= 1:
            mul = mul * (xi - x[i - 1])
        #累加
        res = res + f[i] * mul
    #返回计算结果
    return res
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
#算法5：变步长梯形求积法计算积分

X_INCREMENT = 1e-9    #极小增量   
def cal_fun_x(fun_str,xi): #计算函数fun在xi处的值
    try:
        res = eval(fun_str,None,{'x':xi}) #计算fun表达式的值
    except ZeroDivisionError:   #有除零异常时继续运算xi附近的值
        return cal_fun_x(fun_str,xi + X_INCREMENT)
    #没有异常直接返回计算结果
    return res
    
def var_step_trapezoid(fun_str,a,b,accuracy):   #参数依次为:fun函数（字符串形式），区间[a,b]，精度accuracy
    n = 1	#n为区间被分割的部分数
    h = (b - a) / n  #步长
    t_pre = h / 2 * (cal_fun_x(fun_str,a) + cal_fun_x(fun_str,b))   #t_pre为Tn

    while True:
        t_cur = t_pre / 2   #t_cur为T2n
        for i in range(n):  #h为计算Tn时的步长
            t_cur = t_cur + h / 2 * cal_fun_x(fun_str,a + (i + 0.5) * h)

        if(fabs(t_cur - t_pre) <= accuracy):  #达到给定精度
            return t_cur    #返回当前计算值T2n
        else:
            t_pre = t_cur   #继续迭代
            n = n * 2
            h = (b - a) / n
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
#算法6：龙贝格算法计算积分
    
def romberg(fun_str,a,b,accuracy):  #参数依次为:fun函数字符串，区间[a,b]，精度accuracy
    t_pre = [0]    #t_pre为上一行递推结果
    t_cur = [0,0] #t_cur为当前行递推结果
    k = 0   #k控制递推次数
    n = 1	#n为区间被分割的部分数
    h = (b - a) / n  #h为步长
    #计算T0(0)
    t_pre[0] = h / 2 * (cal_fun_x(fun_str,a) + cal_fun_x(fun_str,b))    
    while True:
        #计算T0(k)
        t_cur[0] = t_pre[0] / 2
        for i in range(n):  #h为计算t_pre的步长
            t_cur[0] = t_cur[0] + h / 2 * cal_fun_x(fun_str, a + (i + 0.5) * h)
        #变化k,n,h
        k = k + 1
        n = n * 2
        h = (b - a) / n
        #递推该行
        for i in range(1,k + 1):
            t_cur[i] = (pow(4,i) * t_cur[i - 1] - t_pre[i - 1]) / (pow(4,i) - 1)
            if(fabs(t_cur[i] - t_cur[i - 1]) <= accuracy):  #精度达到要求
                return t_cur[i]     #返回当前计算结果
        #t_cur赋到t_pre
        t_pre.append(0)
        for i in range(k + 1):
            t_pre[i] = t_cur[i]
        t_cur.append(0)
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
#算法7：改进欧拉法求解一阶常微分方程的初值问题

XY_INCREMENT = 1e-9    #极小增量
def cal_fun_xy(fun_str,xi,yi): #计算函数fun在xi,yi处的值
    try:
        res = eval(fun_str,None,{'x':xi,'y':yi}) #计算fun表达式的值
    except ZeroDivisionError:   #有除零异常时继续运算xi附近的值
        return cal_fun_xy(fun_str,xi + XY_INCREMENT,yi + XY_INCREMENT)
    #没有异常直接返回计算结果
    return res

def improved_euler(fun,a,b,y0,h):   #参数依次为：微分函数f(x,y)的字符串形式,区间[a,b]，初值y0,步长h
    #res数组，第1,2行分别为x,y
    res = [[a],[y0]]
    i = 1	#循环变量
    xi = a + i * h #xi为横坐标值
    while xi <= b:
        #计算xn,yn
        xn = xi - h
        yn = res[1][i - 1]
        y = yn + h * cal_fun_xy(fun,xn,yn) #欧拉法计算yn+1
        y = yn + h / 2 * (cal_fun_xy(fun,xn,yn) + cal_fun_xy(fun,xi,y))   #代入计算新的yn+1
        #结果存至res数组
        res[0].append(xi)
        res[1].append(y)
        i = i + 1
        xi = a + i * h
    #返回二维res数组
    return res
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
#算法8：四阶龙格-库塔求解一阶常微分方程的初值问题

def four_runge_kutta(fun,a,b,y0,h):   #参数依次为：微分函数f(x,y)的字符串形式,区间[a,b]，初值y0,步长h
    #res数组，第1,2行分别为x,y
    res = [[a],[y0]]
    i = 1	#循环变量
    xi = a + i * h #xi为横坐标值
    while xi <= b:
        #计算xn,yn,k1,k2,k3,k4
        xn = xi - h
        yn = res[1][i - 1]
        k1 = cal_fun_xy(fun,xn,yn)
        k2 = cal_fun_xy(fun,xn + h / 2,yn + h / 2 * k1)
        k3 = cal_fun_xy(fun,xn + h / 2,yn + h / 2 * k2)
        k4 = cal_fun_xy(fun,xn + h,yn + h * k3)
        y = yn + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)    #代入公式计算yn+1
        #结果存至res数组
        res[0].append(xi)
        res[1].append(y)
        i = i + 1
        xi = a + i * h
    #返回二维res数组
    return res
#------------------------------------------------------------------------------