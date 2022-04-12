/* lab3.c
 created by Bai
 */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

//global parameteres
int MAX_NUM_OF_ITERS = 100; // 最大运行次数
int MAX_TIME = 30; //max 30 sec comp time. 最大运行时间
int NB_METHOD=1; //0: first descent 1: best descent 方法选择器，0代表first descent, 1代表best descent

/*The following is problem data. For convenience, we
include in the program. For real applications, the data
should be loaded from file(s).*/
int NUM_OF_ITEMS=100;  // 样本数量
int CAPACITY=1000;  // 背包容量

// 各个样本价值
int p[] = {11, 11, 10, 12, 3, 5, 13, 6, 10, 12, 19, 8, 3, 19, 3, 7, 8, 17, 8, 3, 5, 3, 5, 5, 13, 9, 9, 4, 11, 20, 18, 5, 17, 19, 14, 4, 14, 7, 17, 6, 9, 7, 12, 20, 3, 9, 12, 4, 7, 13, 20, 3, 18, 4, 12, 10, 3, 20, 4, 4, 13, 17, 16, 4, 5, 16, 18, 8, 18, 20, 9, 10, 10, 16, 18, 20, 3, 5, 17, 16, 13, 11, 3, 10, 20, 12, 6, 6, 11, 18, 18, 14, 19, 6, 18, 9, 14, 12, 10, 19};
// 各个样本占用的容量
int v[] = {93, 108, 82, 67, 82, 88, 77, 60, 138, 135, 93, 92, 128, 134, 144, 75, 148, 133, 150, 103, 79, 82, 106, 143, 99, 101, 80, 90, 84, 62, 90, 95, 119, 75, 101, 93, 85, 118, 81, 64, 87, 113, 117, 142, 88, 134, 134, 100, 112, 103, 90, 149, 122, 145, 92, 63, 103, 140, 62, 138, 60, 103, 145, 143, 109, 117, 63, 65, 124, 78, 121, 114, 96, 99, 110, 125, 95, 92, 70, 70, 149, 115, 148, 85, 81, 62, 99, 124, 91, 101, 118, 130, 140, 86, 97, 138, 89, 63, 122, 112 };


// 对象: 样本
struct item{
    int indx; //remember the index of items after sorting 样本本来的参数，用来调用数组的值
    int p; //对应样本的价值
    int v; //对应样本的容量
    double ratio;  //对应样本单位容量的价值
};

struct problem{
    int C; //knapsack capacity 背包的容量，即容量上限
    int num_of_items; // 包含样本的数量
    struct item* items; // 背包本身，存放多个样本对象的数组 指针
};

struct solution{
    int objective;  //
    int cap_left;  // 剩余的可分配空间
    struct problem* prob; // 需要解决的问题
    int* x;  // 一个int数组，存放每个样本的使用状况
};

//initialise problem struct with global variables v and p. 初始化问题，以及问题对应的数据样本集
struct problem* init_prob(){
    struct item* my_item = malloc(sizeof(struct item)*NUM_OF_ITEMS);
    for(int i=0; i<NUM_OF_ITEMS; i++){  // 初始化数据集
        my_item[i].indx=i; 
        my_item[i].p=p[i];
        my_item[i].v=v[i];
    }
    struct problem* my_prob = malloc(sizeof(struct problem)); // 初始化问题
    my_prob->C =CAPACITY;
    my_prob->num_of_items = NUM_OF_ITEMS;
    my_prob->items = my_item;
    return my_prob;
}

// 问题解决时，释放内存
void free_problem(struct problem* my_prob)
{
    if(my_prob!=NULL)
    {
        free(my_prob->items); //free item array
        free(my_prob);
    }
}

//create an empty solution 创建空的解决方案
struct solution* create_empty_sol(struct problem* my_prob){
    struct solution* my_sln = malloc(sizeof(struct solution)); 
    my_sln->prob = my_prob;  // 设置对应问题
    my_sln->objective=0; my_sln->cap_left=my_prob->C;  // 设置剩余空间,初始状态为整个背包的容量
    my_sln->x = malloc(sizeof(int)*my_prob->num_of_items);  // 为每个样本的使用状态分配空间
    for(int i=0; i<my_prob->num_of_items; i++) my_sln->x[i]=0; // 初始状态下，所有样本均未被使用
    return my_sln;
}

// 复制解决方案
struct solution* copy_from(struct solution* sln){
    struct solution* my_sln = malloc(sizeof(struct solution));
    my_sln->prob = sln->prob;  //shalow copy 浅复制
    my_sln->objective=sln->objective; my_sln->cap_left=sln->cap_left;
    my_sln->x = malloc(sizeof(int)*sln->prob->num_of_items);
    for(int i=0;i<sln->prob->num_of_items; i++)
        my_sln->x[i] = sln->x[i];
    return my_sln; 
}

// 问题解决时，释放内存
void free_solution(struct solution* sln){
    if(sln!=NULL){
        free(sln->x);
        free(sln);
    }
}

// 打印最佳解决方案，附带正确性检查
void print_sol(struct solution* sln, char* sol_name){
    printf("%s's objective = %d, cap left = %d\n",sol_name, sln->objective, sln->cap_left);
    
    //checking anormaly in the solution struct.
    int cap=0, obj=0;
    for(int i=0; i<sln->prob->num_of_items; i++){
        if(sln->x[sln->prob->items[i].indx]>0){
            obj += sln->prob->items[i].p;
            cap +=sln->prob->items[i].v;
        }
    }
    if(obj != sln->objective) {
        printf("Inconsistent objective value. Correct val = %d\n", obj);
    }
    if(sln->prob->C - cap != sln->cap_left){
        printf("Inconsistent residual capacity value. Correct val =%d\n",
               sln->prob->C- cap);
    }
        
}

// 对比两个样本的单位容量价值
int cmpfunc (const void * a, const void * b) {
    const struct item* item1 = a;
    const struct item* item2 = b;
    if(item1->ratio>item2->ratio) return -1;
    if(item1->ratio<item2->ratio) return 1;
    return 0;
}

// 对比两个样本的序号
int cmpfunc2 (const void * a, const void * b) {
    const struct item* item1 = a;
    const struct item* item2 = b;
    if(item1->indx>item2->indx) return 1;
    if(item1->indx<item2->indx) return -1;
    return 0;
}

//greedy heuristic by packing highest pi/vi 使用样本单位容量价值作为衡量标准的贪心搜索
struct solution* greedy_heuristic(struct problem* prob){
    // 基本异常处理
    if(prob==NULL){
        printf("Missing problem data, please check!\n");
        exit(1);
    }
    // 初始化解决方案
    struct solution* greedy_sln =create_empty_sol(prob);
    int n = prob->num_of_items;
    for(int i=0;i<n; i++){
        // 计算样本单位容量价值
        prob->items[i].ratio = (double)prob->items[i].p/prob->items[i].v;
        //printf("ratio[%d]=%.3f\t",prob->items[i].indx,prob->items[i].ratio);
    }
    //sort items according to ratio 根据样本单位容量价值进行排序
    qsort(prob->items, n, sizeof(struct item), cmpfunc);
    // 初始化问题
    greedy_sln->cap_left=prob->C; greedy_sln->objective=0;
    for(int i=0;i<n; i++){
        //printf("ratio[%d]=%.3f\t",prob->items[i].indx,prob->items[i].ratio);
        // 判断是否还能插入更多对象，如果可以及插入，不行就尝试更小的对象
        if(greedy_sln->cap_left>= prob->items[i].v){
            greedy_sln->x[prob->items[i].indx]=1;
            greedy_sln->cap_left -= prob->items[i].v;
            greedy_sln->objective += prob->items[i].p;
        }
    }
    //make sure original order to items is restored to match index in x 根据样本的原始序号构建解决方案数组
    qsort(prob->items, n, sizeof(struct item), cmpfunc2);
    //print_sol(greedy_sln, "Greedy Solution");
    return greedy_sln;
}

// 找到第一个合适的互换对象并更新最优解
struct solution* first_descent(struct solution* curt_sln)
{//using pair-swaps
    // 存放两个对象的序号
    int item1, item2;
    // 复制现有的解决方案
    struct solution* new_sln = copy_from(curt_sln);
    for(int i=0; i<NUM_OF_ITEMS; i++)
    {
        // 若一个对象被选取了
        if(curt_sln->x[i]>0){
          item1 =i;
          for(int j=0; j<NUM_OF_ITEMS; j++){
            // 取得该对象与另一随机对象的容量
           int v1 = curt_sln->prob->items[item1].v;
              int v2 = curt_sln->prob->items[j].v;
            // 如果该对象不同于随机对象，且该随机对象未被选取，替换后不会导致溢出
          if(i!=j && curt_sln->x[j]==0 &&
             curt_sln->cap_left + v1 >= v2)
          {
              // 将该随机对象的序号进行保存
              item2=j;
              // 计算两个对象的价值差
              int delta =curt_sln->prob->items[item2].p -
                curt_sln->prob->items[item1].p;
              // 如果价值差大于0，将这两个对象对换并将其作为新的解
              if(delta >0){
                  new_sln->x[item1]= 0;
                  new_sln->x[item2]= 1; //swap
                  new_sln->objective=curt_sln->objective + delta;  //delta evaluation 更新最终结果(背包里存放物品的总价值)
                  new_sln->cap_left = new_sln->cap_left - v2 + v1;  // 同时更新背包的剩余容量
                  //print_sol(new_sln, "First Descent Solution");
                  return new_sln;
              }
          }//endif
          }//endfor
        }//endif
    }//endfor
    // 若没有找到合适的替换对象，则返回NULL
    return NULL;
}

// 找到最佳的更换对象并更新最优解
struct solution* best_descent(struct solution* curt_sln)
{//using pair-swaps
    // 存放两个对象的序号
    int item1, item2;
    // 存放最佳的对换对象序号以及其最佳的价值差
    int best_delta = -1000, b_item1=-1, b_item2=-1; //store the best move
    // 复制现有的最佳解
    struct solution* new_sln = copy_from(curt_sln);
    for(int i=0; i<NUM_OF_ITEMS; i++)
    {
        // 若一个对象被选取了
        if(curt_sln->x[i]>0){
          item1 =i;
          // 取得该对象与另一随机对象的容量
          for(int j=0; j<NUM_OF_ITEMS; j++){
              int v1 = curt_sln->prob->items[item1].v;
              int v2 = curt_sln->prob->items[j].v;
              // 如果该对象不同于随机对象，且该随机对象未被选取，替换后不会导致溢出
              if(i!=j && curt_sln->x[j]==0 &&
                 curt_sln->cap_left + v1 >= v2)
              {
                  // 将该随机对象的序号进行保存
                  item2=j;
                  // 计算两个对象的价值差
                  int delta =curt_sln->prob->items[item2].p -
                    curt_sln->prob->items[item1].p;
                     // 如果价值差大于0，将这两个对象对换并将其作为新的解
                  if(delta >0 && delta>best_delta){
                      // 存放当前最佳的互换对象的序号与价值差
                      b_item1 = item1; b_item2 = item2;
                      best_delta= delta;
                  }
              }
          }//endfor
        }//endif
    }//endfor
    // 若存在最佳的互换对象，则将该解决方案作为新的最佳解决方案
    if(best_delta>0){
        new_sln->x[b_item1]= 0;
        new_sln->x[b_item2]= 1; //swap
        new_sln->objective=curt_sln->objective + best_delta;  //delta evaluation
        int v1 = curt_sln->prob->items[b_item1].v;
        int v2 = curt_sln->prob->items[b_item2].v;
        new_sln->cap_left = new_sln->cap_left - v2 + v1;
        //print_sol(new_sln, "Best Descent Solution");
        return new_sln;
    }
    // 若没有，返回NULL
    else
        return NULL;
}
    
// 局部最优算法的主体部分
struct solution* local_search(struct solution* init_sln, int nb_method){
    // 用于存放最佳解法，当前解法
    struct solution *nb_sln, *curt_sln;
    // 初始化当前解决方案
    curt_sln = copy_from(init_sln);
    // 控制运行次数
    for(int k=0; k< MAX_NUM_OF_ITERS; k++)
    {
        // 根据不同的参数，选择不同的解决方案
         if(NB_METHOD==0)
             nb_sln=first_descent(curt_sln);  //first descent local search
         else nb_sln=best_descent(curt_sln);  //best descent local sarch
         if(nb_sln!=NULL)
         {
             // 将最佳解法复制给curt_cln,并释放相应内存
             free_solution(curt_sln);
             curt_sln= copy_from(nb_sln);
             free_solution(nb_sln);
         }
    }
    return curt_sln;
}

// 主函数
int main()
{
    //problem* prob = loaddata("knapsack.txt");
    // 初始化问题
    struct problem* prob = init_prob();
    // 使用贪心算法，找到初始的解决方案
    struct solution* init_sln=greedy_heuristic(prob); //greedy heuristic
    print_sol( init_sln,"Initial Solution");
    
    struct solution *local_optimum_sln;
    // 使用局部搜索算法，找到局部最优解
    local_optimum_sln = local_search(init_sln, NB_METHOD); //local search
    print_sol(local_optimum_sln, "Local Optima Solution");
    
    //free memory for all pointers
    free_problem(prob);
    free_solution(init_sln);
    free_solution(local_optimum_sln);
    return 0;
}

