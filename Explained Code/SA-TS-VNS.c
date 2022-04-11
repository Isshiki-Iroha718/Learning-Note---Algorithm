//
//  main.c
//  mknapsack
//
//  Created by Bai on 14/03/2020.
//  Copyright © 2019 UNNC. All rights reserved.
//

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>
#include <string.h>

/* global parameters */ 
int RAND_SEED[] = {1,20,30,40,50,60,70,80,90,100,110, 120, 130, 140, 150, 160, 170, 180, 190, 200}; // 生成随机数种子
int NUM_OF_RUNS = 2;  //对于每个问题的运行次数
int MAX_TIME = 30;  //max amount of time permited (in sec)   //运行时间限制
int num_of_problems;   //问题个数，这些数据可以从原始文件中读取
enum ALG{SA=1,TS=2,VNS=3};   // 算法选择,1:模拟退火，2:禁忌搜索, 3:变邻域搜索
enum ALG alg=TS;   // 算法选择


/* parameters for evlutionary algorithms*/  // 遗传算法
static int POP_SIZE = 50;   //please modify these parameters according to your problem
int MAX_NUM_OF_GEN = 10000; //max number of generations   总共生成的后代数(代数)
float CROSSOVER_RATE = 0.8;  // 交配概率
float MUTATION_RATE = 0.1;  // 变异概率

/* declare parameters for simulated annealing here */  // 模拟退火算法
float SA_TS =500;    // 初始温度
float SA_TF =10;    // 运行结束时温度标准
float SA_BETA = 0.00000001;   // 改变温度用的参数
int SA_MAX_ITER = 200000000; //total number of runs.  最大运行次数
int SA_ITER_PER_T = 1; //number of runs per temperature 每个温度的运行次数

/* declare parameteres for tabu search here*/ // 禁忌搜索算法
int TABU_TENURE = 20;   // 一个放入禁忌表中的解决方案在表中保留的最大迭代次数
int TS_MAX_ITER = 2000000;    // 算法运行的最大迭代次数
struct tabu_struct{ // 互相交换的两个对象
    int item1;  // 存放两个物品的序号
    int item2;
    int tabu_tenure; //iterations to remain in tabu 一个禁忌解决方案还可以在禁忌表中保留的剩余迭代次数
};

struct move_struct{
    int item1;
    int item2;
};

/* declare parameters for variable neighbourhood search here*/  // 变邻域搜索算法 
int K= 3; // k-opt is used.   邻域的最大数量,也表示最多可以有几个物品参与使用状态交换
int SHAKE_STRENGTH = 12;   // 扰动算子的强度
struct solution_struct best_sln;  //global best solution

//return a random number between 0 and 1  生成随机的0或者1
float rand_01() 
{
    float number;
    number = (float) rand();
    number = number/RAND_MAX;
    //printf("rand01=%f\n", number);
    return number;
}

//return a random nunber ranging from min to max (inclusive) 生成随机整数，范围确定
int rand_int(int min, int max)
{
    int div = max-min+1;
    int val =rand() % div + min;
    //printf("rand_range= %d \n", val);
    return val;
}


struct item_struct{  // 结构体，物品
    int dim; //no. of dimensions   物品的维数
    int* size; //volume of item in all dimensions   每个维度的大小,数据存入一个数组中
    int p;  // 物品的价值
    double ratio;   // 物品对应的性价比
    int indx;  // 物品的序号
};

struct problem_struct{   // 结构体，问题
    int n; //number of items  问题中包含了几个物品
    int dim; //number of dimensions   物品有多少个维度
    struct item_struct* items;  // 存放物品的数组
    int* capacities;  //knapsack capacities  每个背包对应每个维度的剩余容量大小
    double best_obj; //best known objective   当前问题的最佳解决方案得到的最高价值
};

// 当问题被解决时，释放内存
void free_problem(struct problem_struct* prob)
{
    if(prob!=NULL)
    {
        if(prob->capacities !=NULL) free(prob->capacities);
        if(prob->items!=NULL)
        {
            for(int j=0; j<prob->n; j++)
            {
                if(prob->items[j].size != NULL)
                    free(prob->items[j].size);
            }
            free(prob->items);
        }
        free(prob);
    }
}

// 初始化一个问题，传入物品的数量，维数以及接收这个问题的对象my_prob
void init_problem(int n, int dim, struct problem_struct** my_prob)
{
    struct problem_struct* new_prob = malloc(sizeof(struct problem_struct));
    new_prob->n=n; new_prob->dim=dim;
    // 根据物品的数量维度分配内存
    new_prob->items=malloc(sizeof(struct item_struct)*n);
    for(int j=0; j<n; j++)
        new_prob->items[j].size= malloc(sizeof(int)*dim);
    new_prob->capacities = malloc(sizeof(int)*dim);
    *my_prob = new_prob;
}


//example to create problem instances, actual date should come from file 从外部文件中导入问题
struct problem_struct** load_problems(char* data_file)
{
    int i,j,k;
    //int num_of_probs;
    FILE* pfile = fopen(data_file, "r");
    if(pfile==NULL)
        {printf("Data file %s does not exist. Please check!\n", data_file); exit(2); }
    fscanf (pfile, "%d", &num_of_problems);   // 文件读取的第一个数据为问题的数目，数据传入一个全局变量
    
    // 初始化问题，根据全局变量的数据分配相应数目的内存
    struct problem_struct** my_problems = malloc(sizeof(struct problem_struct*)*num_of_problems);
    for(k=0; k<num_of_problems; k++)
    {
        int n, dim, obj_opt;
        fscanf (pfile, "%d", &n); // 传入物品的数量
        fscanf (pfile, "%d", &dim); // 传入物品的维度数量
        fscanf (pfile, "%d", &obj_opt);  // 传入目前的最佳解(背包内存放物品的最大价值)
        
        init_problem(n, dim, &my_problems[k]);  //allocate data memory  
        my_problems[k]->best_obj = obj_opt;  // 初始化最佳选择
        for(j=0; j<n; j++)
        {
            my_problems[k]->items[j].dim=dim;   // 初始化每个物品的维度
            my_problems[k]->items[j].indx=j;   // 初始化每个物品的序号
            fscanf(pfile, "%d", &(my_problems[k]->items[j].p)); //read profit data  读取每个物品的价值
            //printf("item[j].p=%d\n",my_problems[k]->items[j].p);
        }
        // 读取每个物品在各个维度上的大小
        for(i=0; i<dim; i++)
        { 
            for(j=0; j<n; j++)  
            {
                fscanf(pfile, "%d", &(my_problems[k]->items[j].size[i])); //read size data  
                //printf("my_problems[%i]->items[%i].size[%i]=%d\n",k,j,i,my_problems[k]->items[j].size[i]);
            }
        }
        // 读取每个维度的限制条件
        for(i=0; i<dim; i++){
            fscanf(pfile, "%d", &(my_problems[k]->capacities[i]));   
            //printf("capacities[i]=%d\n",my_problems[k]->capacities[i] );
        }
    }
    fclose(pfile); //close file
    return my_problems;
}

// 解决方案的结构体
struct solution_struct{
    struct problem_struct* prob; //maintain a shallow copy of the problem data  对一个问题内容的复制
    float objective;  // 背包中存放物品的最大价值
    int feasibility; //indicate the feasiblity of the solution  衡量解法的可行性
    int* x; // solution encoding vector  使用一个数组保存解法中各个物品的使用情况
    int* cap_left; //capacity left in all dimensions  使用一个数组保存每个维度的剩余空间
};

// 对于解决方案，释放内存
void free_solution(struct solution_struct* sln)
{
    if(sln!=NULL)
    {
        free(sln->x);
        free(sln->cap_left);
        sln->objective=0;
        sln->prob=NULL;
        sln->feasibility=false;
    }
}

//copy a solution from another solution  复制对应解法
bool copy_solution(struct solution_struct* dest_sln, struct solution_struct* source_sln)
{
    if(source_sln ==NULL) return false;
    if(dest_sln==NULL)
    {
        dest_sln = malloc(sizeof(struct solution_struct));
    }
    else{
        free(dest_sln->cap_left);
        free(dest_sln->x);
    }
    int n = source_sln->prob->n;
    int m =source_sln->prob->dim;
    dest_sln->x = malloc(sizeof(int)*n);
    dest_sln->cap_left=malloc(sizeof(int)*m);
    for(int i=0; i<m; i++)
        dest_sln->cap_left[i]= source_sln->cap_left[i];
    for(int j=0; j<n; j++)
        dest_sln->x[j] = source_sln->x[j];
    dest_sln->prob= source_sln->prob;
    dest_sln->feasibility=source_sln->feasibility;
    dest_sln->objective=source_sln->objective;
    return true;
}


// 评估解法的最优性
void evaluate_solution(struct solution_struct* sln)
{
    //evaluate the feasibility and objective of the solution  
    sln->objective =0; sln->feasibility = 1;  //
    struct item_struct* items_p = sln->prob->items;
    
    for(int i=0; i< items_p->dim; i++)
    {
        sln->cap_left[i]=sln->prob->capacities[i];
        for(int j=0; j<sln->prob->n; j++)
        {
            sln->cap_left[i] -= items_p[j].size[i]*sln->x[j];
            if(sln->cap_left[i]<0) {
                sln->feasibility = -1*i; //exceeding capacity
                return;
            }
        }
    }
    if(sln->feasibility>0)
    {
        for(int j=0; j<sln->prob->n; j++)
        {
            sln->objective += sln->x[j] * items_p[j].p;
        }
    }
}

//output a given solution to a file 将输出写入一个文件中
void output_solution(struct solution_struct* sln, char* out_file)
{
    if(out_file !=NULL){
        FILE* pfile = fopen(out_file, "a"); //append solution data
        double gap=1000; //set to an arbitrarily large number if best known solution is not availabe.
        if(best_sln.prob->best_obj!=0) gap=  100*(best_sln.prob->best_obj - best_sln.objective)/best_sln.prob->best_obj;
        fprintf(pfile, "%i \t %0.2f\n", (int)sln->objective, gap);
        for(int i=0; i<sln->prob->n; i++)
        {
            fprintf(pfile, "%i ", sln->x[i]);
        }
        fprintf(pfile, "\n");
        /*for(int j=0; j<sln->prob->n; j++)
            fprintf(pfile, "%i ", sln->prob->items[j].p);
        fprintf(pfile, "\n");*/
        fclose(pfile);
    }
    else
        printf("sln.feas=%d, sln.obj=%f\n", sln->feasibility, sln->objective);
}


//check the feasiblity and obj values of solutons from solution_file.
//return 0 is all correct or output infeasible solutions to files.
int check_solutions(struct problem_struct** my_problems, char* solution_file)
{
    FILE * pfile= fopen(solution_file, "r");
    if(pfile==NULL)
    {
        printf("Solution file %s does not exist. Please check!\n", solution_file);
        exit(2);
    }
    float val_obj;
    int val;
    fscanf (pfile, "%i", &val);
    if(val != num_of_problems)
    {
        printf("The stated number of solutions does not match the number of problems.\n");
        exit(3);
    }
    
    fscanf(pfile, "%d", &val); //get rid of gap information.
    struct solution_struct temp_sln;
    int count=0, k=0;
    int n, dim;
    while(fscanf (pfile, "%f", &val_obj)!=EOF && k<num_of_problems)
    {
        //val_obj = val;
        n= my_problems[k]->n;  dim= my_problems[k]->dim;
        temp_sln.x = malloc(sizeof(int)*n);
        temp_sln.cap_left=malloc(sizeof(int)*dim);
        temp_sln.prob = my_problems[k];
        while(fscanf (pfile, "%i", &val)!=EOF)
        {
            if(val<0 || val>1) {fclose(pfile);  return k+1;} //illeagal values
            temp_sln.x[count] = val;
            count++;
            if(count==n)
            {
                evaluate_solution(&temp_sln);
                if(!temp_sln.feasibility || fabs(temp_sln.objective - val_obj)>0.01)
                {
                    fclose(pfile);
                    //printf("feasb=%i, obj= %f, val=%i\n",temp_sln.feasibility, temp_sln.objective, val_obj);
                    //output_solution(&temp_sln, "my_debug.txt");
                    return k+1;  //infeasible soltuion or wrong obj
                }
                else{
                    break;
                }
            }
        }
        count=0; k++;
        
        free(temp_sln.x); free(temp_sln.cap_left);
    }
    fclose(pfile);
    return 0;
}

//modify the solutions that violate the capacity constraints
void feasibility_repair(struct solution_struct* pop)
{
    //todo
}

//local search
void local_search_first_descent(struct solution_struct* pop)
{
    //todo
}


//update global best solution from sln  更新最佳的解法
void update_best_solution(struct solution_struct* sln)
{
    if(best_sln.objective < sln->objective)
    copy_solution(&best_sln, sln);
}

// 比较两个物品的性价比
int cmpfunc1(const void* a, const void* b){
    const struct item_struct* item1 = a;
    const struct item_struct* item2 = b;
    if(item1->ratio>item2->ratio) return -1;
    if(item1->ratio<item2->ratio) return 1;
    return 0;
    }

// 比较两个物品的序号
int cmpfunc2 (const void * a, const void * b) {
        const struct item_struct* item1 = a;
        const struct item_struct* item2 = b;
        if(item1->indx>item2->indx) return 1;
        if(item1->indx<item2->indx) return -1;
        return 0;
    }

// 比较两个解法的价值
int cmpfunc_sln (const void * a, const void * b) {
    const struct solution_struct* sln1 = a;
    const struct solution_struct* sln2 = b;
    if(sln1->objective > sln2 ->objective) return -1;
    if(sln1->objective < sln2 ->objective) return 1;
    return 0;
}

//a greedy heuristic solution based on profit/volume ratio  基于性价比，使用贪心算法得出当前最佳的解法
struct solution_struct* greedy_heuristic(struct problem_struct* prob)
{
    for(int i=0; i<prob->n;i++){
        // 针对每个物品,计算性价比的方式
        double avg_size=0;
        struct item_struct* item_i = &prob->items[i];
        for(int d=0; d< prob->dim; d++){
            avg_size += (double)item_i->size[d]/prob->capacities[d];
        }
        item_i->ratio = item_i->p/avg_size;
    }
    qsort(prob->items, prob->n, sizeof(struct item_struct), cmpfunc1);
    
    // 初始化解法
    struct solution_struct* init_sln = malloc(sizeof(struct solution_struct));
    init_sln->prob=prob;    init_sln->objective =0;
    init_sln->x = malloc(sizeof(int)*prob->n);
    init_sln->cap_left = malloc(sizeof(int)*prob->dim);
    int* cap = malloc(sizeof(int)*prob->dim);
    int i=0, d=0;
    for(d=0; d<prob->dim; d++) cap[d]=0; //aggregated volume
    for(i=0; i<prob->n; i++)
    {
        struct item_struct* item_i = &prob->items[i];
        //printf("item[%d].ratio = %.3f\t",item_i->indx,prob->items[i].ratio);
        for(d=0; d<prob->dim; d++){
            if(cap[d] + item_i->size[d] > prob->capacities[d])
                break; //infeasible to pack this item, try next
        }
        if(d>=prob->dim){
            init_sln->x[item_i->indx] = 1;
            init_sln->objective += item_i->p;
            for(d=0; d<prob->dim; d++){
                cap[d] += item_i->size[d];
            }
            //printf("packing item %d\n", item_i->indx);
        }
        else init_sln->x[item_i->indx] =0;
    }
    for(d=0; d<prob->dim; d++){
        init_sln->cap_left[d] = prob->capacities[d]- cap[d];
    }
    free(cap);
    //restore item original order by sorting by index.
    qsort(prob->items, prob->n, sizeof(struct item_struct), cmpfunc2);
    
    evaluate_solution(init_sln);
    //output_solution(init_sln, "greedy_sln.txt");
    printf("Init_sln obj=\t%0.0f\tfeasiblity = %d.\n", init_sln->objective, init_sln->feasibility);
    return init_sln;
}

//Simulated Annealing 模拟退火算法的实现
int SimulatedAnnealing(struct problem_struct* prob)
{// 计时用的参数
    clock_t time_start, time_fin;
    time_start = clock();
    double time_spent=0;
    int iter =0;    // 迭代次数
    double temperature=SA_TS;   // 初始化温度
    best_sln.prob = prob;   // 初始化问题的解决方法
    struct solution_struct* curt_sln = greedy_heuristic(prob);  // 最初的解法是根据贪心算法得到的解法
    update_best_solution(curt_sln);     // 将目前的解法作为最佳解法
    struct solution_struct* rand_neighb=malloc(sizeof(struct solution_struct));  // 初始化相邻解决方案，随机生成
    rand_neighb->cap_left= malloc(sizeof(int)*prob->dim);   // 根据维数和每个维度的容量分配内存
    rand_neighb->x = malloc(sizeof(int)*prob->n);   // 根据物品的数量分配内存
    while(iter<SA_MAX_ITER && time_spent < MAX_TIME && temperature > SA_TF)   // 循环停止条件:达到最大迭代次数，达到最大运行时间，达到退火最低温度
    {
        //add your SA code here
        copy_solution(rand_neighb, curt_sln);
        int item1, item2;
        // item1 存放一个被解决方案选择的物品, item2 存放一个未被选择的物品
        item1 = rand_int(0, prob->n-1);   // 生成随机数(0,物品的数量-1) 实际上是物品的参数
        if(curt_sln->x[item1] ==1){         // 如果该物品确定在这个解法中被选中了
            item2 = rand_int(0, prob->n-1);     // 生成随机数(0,物品的数量-1) 实际上是物品的参数
            while(curt_sln->x[item2] ==1){//careful, this might lead to deadloop
                item2 = rand_int(0, prob->n-1);   
            }
        }
        else{  // 若没有被选中
            item2 = rand_int(0, prob->n-1);  // 随机找一个被当前解决方案选择的物品
            while(curt_sln->x[item2] ==0){//careful, this might lead to deadloop
                item2 = rand_int(0, prob->n-1);  
            }
            int temp = item1;  // 交换二者的选择状态,使得item1与item2的定义保持原状
            item1 = item2;
            item2 = temp;
        }
        
        //testing potential constraint violations after swap  检测潜在的溢出现象(交换之后导致部分维度的capacity超出限制条件)
        bool flag=true; // 检测用的标记
        for(int d=0; d<prob->dim; d++){  // 对每个维度都进行一次检测
            if(rand_neighb->cap_left[d] + prob->items[item1].size[d] <
               prob->items[item2].size[d]){
                flag=false;     // 出现异常就将标记改为false
                break;
            }
        }
        if(flag){//can swap  // 如果标记是正常的
            float delta = prob->items[item2].p - prob->items[item1].p;   // 计算交换前后解法得到的解是否是最优的
            // 如果更新后的解决方案更优，无条件接受该方案
            // 如果没有做到更优，则按照一定的概率接受这个解决方案
            if(delta>=0 || (delta<0 && exp(delta/temperature)> rand_01())){ 
                // 接受解法之后更新当前解决方案的相应数据   
                rand_neighb->x[item1]=0;
                rand_neighb->x[item2]=1;
                rand_neighb->objective += delta;
                // 更新每个维度的剩余容量capacity
                for(int d=0; d<prob->dim; d++){
                    rand_neighb->cap_left[d] +=  prob->items[item1].size[d] - prob->items[item2].size[d];
                }
            }
            // 将得到的解决方案复制给curt_sln,并将其作为当前的最优解,虽然其甚至可能不是局部最优解
            copy_solution(curt_sln, rand_neighb);
            update_best_solution(curt_sln);   // 如果该解不如之前的解决方案,这步骤将不进行任何操作(即不会修改最优解，原有的最优解将会保留)
            if(iter%100 ==0)  // 没迭代100次，输出一次当前温度与其他相应数据
                printf("tempereature=%0.2f, curt obj =%0.0f,\t best obj=%0.0f\n",temperature, curt_sln->objective, best_sln.objective);
        }
        temperature = temperature/(1+SA_BETA*temperature);   // 降温步骤
        iter++;     // 记录迭代次数
        time_fin=clock();   // 记录运行终止时间
        time_spent = (double)(time_fin-time_start)/CLOCKS_PER_SEC;   //计算总运行时间
    }

    //output_solution(&best_sln, "SA_debug.txt");
    free_solution(curt_sln);
    free_solution(rand_neighb);
    free(rand_neighb);
    return 0;
}

// 根据解法判断能否将两个物品切换
bool can_swap(struct solution_struct* sln, int out, int in)
{
    for(int d =0; d<sln->prob->dim; d++)
    {
        if(sln->cap_left[d]+sln->prob->items[out].size[d] < sln->prob->items[in].size[d]) // 检查每个维度的容量是否溢出
            return false;
    }
    return true;
}

struct move_struct get_move(struct solution_struct* new_sln, struct solution_struct* orig_sln){
    struct move_struct mv;
    int count=0;  // 存放状态的参数
    for(int i=0; i<orig_sln->prob->n; i++)
    {
        // 如果新解法对于某个物品的使用状况与原有解法不同
        if(new_sln->x[i] != orig_sln->x[i]){
            // 将第一个不同状态的物品作为第一个改变对象
            if(count==0){ mv.item1 = i; count =1;}
            // 并将最后一个不同状态的物品作为第二个改变对象
            else mv.item2=i;
        }
    }
    return mv;
}

//return the position of the move in the tabu lift. return -1 if not in the list  取得一对物品在指定禁忌表中的位置
int get_tabu_pos(struct tabu_struct* tabu_list, struct move_struct mv, int count){
    for(int i=0; i<count; i++)
    {
        if(mv.item1==tabu_list[i].item1 && mv.item2 == tabu_list[i].item2) return i;
        if(mv.item1==tabu_list[i].item2 && mv.item2 == tabu_list[i].item1) return i;
    }
    return -1;
}

// 更新禁忌表,对应的是其中的一对物品
void update_tabu_list(struct tabu_struct* tabu_list, struct move_struct mv, int count){
    // 在进行插入操作前,可以确定表内所有的元素都已经经过一次周期,所以要更新表内所有的元素(剩余周期 - 1)
    //update  如果一个元素的禁忌周期达到阈值，将该元素移出禁忌表
    for(int i=0; i<count; i++)
    {
        // 对于禁忌表中的物品要减去其禁忌周期
        tabu_list[i].tabu_tenure -= 1;
        if(tabu_list[i].tabu_tenure<=0)
        {
            tabu_list[i].item1=-1; tabu_list[i].item2=-1; //remove from tabu
        }
    }
    // 在插入操作之前先进行更新操作,确保禁忌表中过期的元素可以先行移出
    //insert mv 
    int p=get_tabu_pos(tabu_list, mv, count);   // 判断一对物品是否被储存在指定禁忌表内,并得到这对物品位于禁忌表中的位置
    if(p>=0)
    {//若该对物品已经被存放在禁忌表中,则更新其禁忌周期(禁忌期限)至刚置入禁忌表时的状态
        tabu_list[p].tabu_tenure= TABU_TENURE;
    }
    else{//若该对物品之前没有放在禁忌表中,找到禁忌表中的第一个空位并将其置入
        for(int i=0; i<count; i++)
        {
            if(tabu_list[i].item1<0 || tabu_list[i].item2<0) { //insert in first empty entry
                tabu_list[i].item1=mv.item1; tabu_list[i].item2=mv.item2;
                tabu_list[i].tabu_tenure = TABU_TENURE;
                break;
            }
        }
    }
}

//TS method 禁忌算法的实现
int TabuSearch(struct problem_struct* prob){
    clock_t time_start, time_fin;
    time_start = clock();
    double time_spent=0;
    int iter =0;  // 初始化迭代次数
    struct tabu_struct tabu_list[100];  //arbitrially large enough  设置一个大小为100的禁忌表
    for(int i=0; i<100; i++) {  // 初始化禁忌表内容，默认状态为空，即用不可用的数字来表示空
        tabu_list[i].item1=-1; tabu_list[i].item2=-1; tabu_list[i].tabu_tenure=0;
    }
    best_sln.prob = prob;
    struct solution_struct* curt_sln = greedy_heuristic(prob);  // 初始化的解决方案由贪心算法生成
    update_best_solution(curt_sln);  // 将目前的最优解决方案改为初始化方案
    // 构造当前解决方案的邻域,一个领域包含多个解法
    int neighb_size_max = prob->n*prob->n/4; //maximum num of neighbouring solutions  邻域大小最大为(物品数量)²/4
    struct solution_struct* N_S=malloc(sizeof(struct solution_struct)*neighb_size_max);   // 初始化邻域,根据邻域的大小分配内存
    // 初始化解法的各个内容的内存以及初始数据
    for(int k=0; k<neighb_size_max; k++)
    {
        N_S[k].cap_left= malloc(sizeof(int)*prob->dim);
        N_S[k].x = malloc(sizeof(int)*prob->n);
        copy_solution(&N_S[k], curt_sln);
        N_S[k].objective = -10000;
    }
    
    // 算法循环终止条件: 运行时间达到阈值，迭代次数达到阈值
    // 注意,每次迭代时不会将禁忌表置空,但会重置所有的邻域(虽然在每次迭代结束之后没有删除邻域的操作)
    while(iter< TS_MAX_ITER && time_spent < MAX_TIME)
    {
        int count=0; //count the number of feasible neighbours 初始状态下没有邻域
        int i, j;   // 用于表示可交换的两个物品
        // 每次生成邻域使更换不同的两个物品(按照物品对来匹配)
        for(i=0; i<prob->n; i++){
            if(curt_sln->x[i]>0){
                //item1 =i;
                for(j=0; j<prob->n; j++){
                    //item2 = j;
                    // 交换其中两个物品的使用状态并生成相应的邻域
                    // 选出的第二件物品必须是没有被使用的
                    if(i!=j && curt_sln->x[j]==0 && can_swap(curt_sln,i, j))
                    {// 将当前解决方案复制给邻域并进行操作
                        copy_solution(&N_S[count], curt_sln);
                        // 对每一个维度的剩余容量进行操作
                        for(int d=0; d<prob->dim; d++){
                            N_S[count].cap_left[d] = N_S[count].cap_left[d]+ prob->items[i].size[d]-prob->items[j].size[d];
                        }
                        // 与模拟退火算法操作类似，计算差值并进行交换
                        int delta =curt_sln->prob->items[j].p -curt_sln->prob->items[i].p;
                        N_S[count].objective = curt_sln->objective + delta;
                        N_S[count].x[i] = 0; //swap
                        N_S[count].x[j] = 1;
                        /*evaluate_solution(&N_S[count]);
                        if(N_S[count].feasibility<0){
                            printf("infeasible sln=%d\n", N_S[count].feasibility);
                            output_solution(&N_S[count], "TS_debug.txt");
                        }*/
                        // 上述操作结束之后，邻域的总数量增加1个
                        // 再下一次迭代中，可以操作下一个邻域
                        count++;
                    }
                }
            }
        }
        // 对于各个邻域的性价比进行排序
        qsort(N_S, count, sizeof(struct solution_struct), cmpfunc_sln);
        for(int k =0; k<count; k++){ // 对每个邻域进行如下操作
            //printf("sln[%d].obj = %.0f\n",count,N_S[k].objective);
            struct move_struct mv = get_move(&N_S[k], curt_sln);  // 找到每个邻域可对换的一对物品
            if(get_tabu_pos(tabu_list, mv, 100)>=0)  // 此处的100指的是禁忌表的大小,如果该对物品位于禁忌表中
            {
                if(N_S[k].objective > best_sln.objective){ // 如果该解决方案比目前的好
                    copy_solution(curt_sln, &N_S[k]); //aspiration criteria 渴望准则,将禁忌表中的该解法解禁已获得更优的解法
                    update_tabu_list(tabu_list, mv,100);  // 更新禁忌表
                    break;
                }
            }
            else{  //若该物品尚未存在禁忌表中
                copy_solution(curt_sln, &N_S[k]);  // 直接存入禁忌表,并将该解决方案作为当前最佳解
                update_tabu_list(tabu_list, mv,100);
                break;
            }
        }
        update_best_solution(curt_sln); // 检测当前的解决方案是否为最佳的解决方案,如果是则更新,不是则不进行任何操作
        if(iter%100 ==0)
        printf("iter=%d, curt obj =%0.0f,\tbest obj=%0.0f\n",iter, curt_sln->objective, best_sln.objective);
        // 计算运行时间,增加迭代次数
        iter++;
        time_fin=clock();
        time_spent = (double)(time_fin-time_start)/CLOCKS_PER_SEC;
    }
    // 释放内存操作
    output_solution(&best_sln, "TS_debug.txt");
    free_solution(curt_sln);
    for(int k=0; k<neighb_size_max; k++)
        free_solution(&N_S[k]);
    free(N_S);
    return 0;
}

// 检测物品是否可以改变状态
// nb_indx 表示有几个物品参与改变状态
// move 作为一个数组,存放改变状态的物品的序号
bool can_move(int nb_indx, int* move, struct solution_struct* curt_sln ){
    bool ret=true;
    // 交换的物品数量为1
    // 默认将这个物品从未使用状态转化为使用状态
    // 从使用该函数的位置可以看出这一点
    if(nb_indx==1)
    {
        int i = move[0];
        if(i<0) return false; // 说明这个物品不存在
        for(int d=0; d<curt_sln->prob->dim; d++){
            if(curt_sln->cap_left[d] < curt_sln->prob->items[i].size[d]) // 如果当这个物品换进时产生了容量溢出
                return false;
        }
    }
    else if(nb_indx==2){  // 交换的物品数量为2
        ret=can_swap(curt_sln, move[0], move[1]);  // 直接使用现成的函数可以解决这个问题
    }
    else if(nb_indx==3){//3-item swap  三个物品参与交换，这里考虑了两种情况
        int i= move[0], j= move[1], k= move[2];   // i一定表示换出的物品,k一定表示换入的物品
        if(i<0 || j<0 || k<0) return false;   // 检测物品是否存在
        if(curt_sln->x[j]>0) {//2-1 swap  判断j是换入还是换出,此处判断j是换出
            for(int d=0; d<curt_sln->prob->dim; d++){
                if(curt_sln->cap_left[d] + curt_sln->prob->items[i].size[d] +
                   curt_sln->prob->items[j].size[d] < curt_sln->prob->items[k].size[d])  //看麻了,解释一下
                   // 剩余空间+物品i的占用空间+物品j的占用空间 < 物品k的占用空间
                   // 以上占用空间均为维度d的占用空间
                    return false;
            }
        }
        else {//1-2 swap   此处判断j是换入
            for(int d=0; d<curt_sln->prob->dim; d++){
                if(curt_sln->cap_left[d] + curt_sln->prob->items[i].size[d] <
                   curt_sln->prob->items[j].size[d] + curt_sln->prob->items[k].size[d])
                    return false;
            }
        }
    }
    // 传进去的nb_indx比3要大
    // 教授偷懒,不想写超过3个物品的交换状态算法
    // 不过有一说一也没啥必要
    else ret=false;
    return ret;
}

// 根据参数执行相应改变状态操作
// 该函数会修改sln,即得到一个新的解决方案
// 返回true表示操作成功执行
bool apply_move(int nb_indx, int* move, struct solution_struct* sln ){
    bool ret=true;
    if(nb_indx==1)  // 对1个物品进行操作
    {
        int i = move[0];
        if(i<0) return false;
        for(int d=0; d<sln->prob->dim; d++){
            sln->cap_left[d] -= sln->prob->items[i].size[d]; // 换入物品后,剩余空间减少
        }
        sln->objective += sln->prob->items[i].p;  // 解法的价值提升
        sln->x[i]=1;  // 该物品的使用状态设置为1,即正在被使用
        
        //printf("success\n");
    }
    else if(nb_indx==2){  // 对2个物品进行操作
        for(int d=0; d<sln->prob->dim; d++){  
            sln->cap_left[d] = sln->cap_left[d] + sln->prob->items[move[0]].size[d]-
                sln->prob->items[move[1]].size[d];
        }
        sln->objective += sln->prob->items[move[1]].p-sln->prob->items[move[0]].p;
        sln->x[move[0]]=0; sln->x[move[1]]=1;
    }
    else if(nb_indx==3){//3-item swap  对3个物品进行操作
        int i= move[0], j= move[1], k= move[2];
        if(i<0 || j<0 || k<0) return false;
        if(sln->x[j]>0) {//2-1 swap   两种不同的操作手段,根据物品j的使用状况来定
            for(int d=0; d<sln->prob->dim; d++){
                sln->cap_left[d] = sln->cap_left[d]+sln->prob->items[i].size[d] +
                    sln->prob->items[j].size[d] - sln->prob->items[k].size[d];
            }
            sln->objective += sln->prob->items[k].p-sln->prob->items[i].p-sln->prob->items[j].p;
            sln->x[i]=0; sln->x[j]=0; sln->x[k]=1;
        }
        else {//1-2 swap
            for(int d=0; d<sln->prob->dim; d++){
                sln->cap_left[d] = sln->cap_left[d]+sln->prob->items[i].size[d] -
                    sln->prob->items[j].size[d] - sln->prob->items[k].size[d];
            }
            sln->objective += sln->prob->items[j].p+sln->prob->items[k].p-sln->prob->items[i].p;
            sln->x[i]=0; sln->x[j]=1; sln->x[k]=1;
        }
        
    }
    // 上述操作与can_move极其类似,这里就不再赘述了
    else ret=false;
    return ret;
}

//nb_indx <=3
struct solution_struct* best_descent_vns(int nb_indx, struct solution_struct* curt_sln)
{ // 用于寻找最佳邻域
    struct solution_struct* best_neighb = malloc(sizeof(struct solution_struct));
    best_neighb->cap_left = malloc(sizeof(int)*curt_sln->prob->dim);
    best_neighb->x = malloc(sizeof(int)*curt_sln->prob->n);
    copy_solution(best_neighb, curt_sln);
    int n=curt_sln->prob->n;   // 问题中包含物品的数量
    int curt_move[] ={-1,-1,-1}, best_move []={-1,-1,-1}, delta=0, best_delta=0;  //storing best neighbourhood moves 记录邻域移动
    switch (nb_indx)
    {
        case 1: //check whether any items can be inserted.
            for(int i=0; i<n; i++){
                if(curt_sln->x[i]>0) continue;  // 不可以把已经在使用的物品改为使用状态
                curt_move[0]=i;  // 向数组中填充移动物品的序号
                if(can_move(nb_indx, &curt_move[0], best_neighb)){   // 如果可以移动,则计算移动后与当前最佳解决方案价值的差值
                    delta = curt_sln->prob->items[i].p;
                    if(delta> best_delta) { // 如果此时的差值最大,则将其作为当前最佳
                        best_delta = delta; best_move[0] = i;
                    }
                }
            }
            // 如果该邻域存在更好的解决方案,则执行这次移动
            if(best_delta>0) {    apply_move(nb_indx, &best_move[0], best_neighb);}
            break;
        case 2:
            for(int i=0; i<n; i++){ // 与Case 1一样的操作流程
                if(curt_sln->x[i]<=0) continue;
                for(int j=0; j<n; j++){
                    if(curt_sln->x[j]==0)
                    {
                        curt_move[0]= i; curt_move[1]= j; curt_move[2]=-1; // 不存在第三个参与交换的物品
                        if(can_move(nb_indx, &curt_move[0], best_neighb)){
                            delta = curt_sln->prob->items[j].p -curt_sln->prob->items[i].p;
                            if(delta > best_delta){
                                best_delta = delta; best_move[0] = i; best_move[1] = j; best_move[2]=-1;
                            }
                        }
                    }
                }
            }
            if(best_delta>0) { apply_move(nb_indx, &best_move[0], best_neighb);}
            break;
        case 3:  // 三个物品交换状态时,考虑两种交换情况
            //2-1 swap
            for(int i=0; i<n; i++){
                if(curt_sln->x[i]==0) continue;
                for(int j=0; j!=i&&j<n; j++){
                    if(curt_sln->x[j]==0) continue;
                    for(int k=0;k<n;k++){
                        if(curt_sln->x[k] == 0)
                        {
                            curt_move[0]=i; curt_move[1]=j; curt_move[2]=k;
                            if(can_move(nb_indx, &curt_move[0], best_neighb)){
                                delta = curt_sln->prob->items[k].p -curt_sln->prob->items[i].p-curt_sln->prob->items[j].p;
                                if(delta > best_delta){
                                    best_delta = delta; best_move[0] = i; best_move[1] = j; best_move[2]=k;
                                }
                            }
                        }
                    }
                }
            }
            //1-2 swap
            for(int i=0; i<n; i++){
                if(curt_sln->x[i]==0) continue;
                for(int j=0; j<n; j++){
                    if(curt_sln->x[j]>0) continue;
                    for(int k=0;k!=j&&k<n;k++){
                        if(curt_sln->x[k] == 0)
                        {
                            curt_move[0]=i; curt_move[1]=j; curt_move[2]=k;
                            if(can_move(nb_indx, &curt_move[0], curt_sln)){
                                delta = curt_sln->prob->items[k].p +curt_sln->prob->items[j].p-curt_sln->prob->items[i].p;
                                if(delta > best_delta){
                                    best_delta = delta; best_move[0] = i; best_move[1] = j; best_move[2]=k;
                                }
                            }
                        }
                    }
                }
            }
            if(best_delta>0) { apply_move(nb_indx, &best_move[0], best_neighb);}
            break;
        default:  // 不支持过多的物品交换状态操作,最多3个
            printf("Neighbourhood index is out of the bounds, nothing is done!\n");
    }
    return best_neighb;
}

void vns_shaking(struct solution_struct* sln, int strength)
{//using random pair-wise swap  强度最多到12,根据设置在文件头部的扰动强度的常数
    int n= sln->prob->n;   // 问题包含物品的数量
    int m =0, try=0;   // 尝试次数最多200回,扰动强度不能大于设定的strength
    while(m<strength && try<200)
    {
        // 随机将两个物品的使用状况
        int move[2];
        int i = rand_int(0, n-1);
        if(sln->x[i]>0){  // 如果随机选出的一个物品的使用状况是在被使用
            move[0] = i;
            do{
                move[1] =rand_int(0, n-1);  // 找到另一个没有使用的物品
            }while(sln->x[move[1]]>0);
        }
        else{  // 如果随机选出的一个物品的使用状况是没有在被使用
            move[1]=i;   
            do{
                move[0] = rand_int(0, n-1);  // 找到另一个在使用的物品
            }while(sln->x[move[0]]<=0);
        }
        if(can_swap(sln, move[0], move[1]))   // 如果两个物品可以交换状态
        {
            apply_move(2, &move[0], sln);   // 对其实行交换并提高扰动强度
            m++; 
        }
        try++;  // 一次尝试的完成,计数用参数
    }
}
//VNS 变邻域搜索算法
int varaible_neighbourhood_search(struct problem_struct* prob){
    clock_t time_start, time_fin;
    time_start = clock();
    double time_spent=0;
    int nb_indx =0; //neighbourhood index 邻域参数,也可以用来表示改变状态的物品数目(含义就是这个)
    
    best_sln.prob = prob;
    // 生成初始解,并将该解作为当前的最佳解
    struct solution_struct* curt_sln = greedy_heuristic(prob);
    update_best_solution(curt_sln);
    
    int shaking_count =0; // 计算扰动的次数,除了计数外没发现有什么用
    while(time_spent < MAX_TIME) //note that final computational time can be way beyond the MAX_TIME if best_descent is time consuming VNS算法唯一的终止条件就是运行时间
    {
        while(nb_indx<K){  // nb_indx<3
            struct solution_struct* neighb_s=best_descent_vns(nb_indx+1, curt_sln); //best solution in neighbourhood nb_indx
            if(neighb_s->objective > curt_sln->objective){ // 当邻域解法比现有解法要好时
                copy_solution(curt_sln, neighb_s); // 将该解法作为现在的最优解法
                nb_indx=1;
            }
            else nb_indx++;  // 当在领域中的解法均没有比现有解法要好,就去下一个邻域搜索
            free_solution(neighb_s);free(neighb_s);
        }
        update_best_solution(curt_sln);  // 检测得出的现有解法是否可以作为最优解法
        double gap=1000; //set to an arbitrarily large number if best known solution is not availabe. 如果没有得到最高价值的解法,设置1000为默认值
        if(best_sln.prob->best_obj!=0) gap=  100*(best_sln.prob->best_obj - best_sln.objective)/best_sln.prob->best_obj;  // 如果最高价值的解法,将其作为gap的值
        printf("shaking_count=%d, curt obj =%0.0f,\t best obj=%0.0f,\t gap= %0.2f%%\n",shaking_count, curt_sln->objective, best_sln.objective, gap);
        vns_shaking(curt_sln, SHAKE_STRENGTH); //shaking at a given strength. This can be made adaptive
        //vns_shaking(curt_sln, shaking_count/100+1); //re-active shaking
        shaking_count++;
        nb_indx=0;
        
        time_fin=clock();
        time_spent = (double)(time_fin-time_start)/CLOCKS_PER_SEC;
    }

    //output_solution(&best_sln, "vns_results.txt");
    free_solution(curt_sln); free(curt_sln);
    return 0;
}

// 主要操作函数,懒得管了就这样吧
int main(int argc, const char * argv[]) {
    
    printf("Starting the run...\n");
    char data_file[50]={"somefile"}, out_file[50], solution_file[50];  //max 50 problem instances per run
    if(argc<3)
    {
        printf("Insufficient arguments. Please use the following options:\n   -s data_file (compulsory)\n   -o out_file (default my_solutions.txt)\n   -c solution_file_to_check\n   -t max_time (in sec)\n");
        return 1;
    }
    else if(argc>9)
    {
        printf("Too many arguments.\n");
        return 2;
    }
    else
    {
        for(int i=1; i<argc; i=i+2)
        {
            if(strcmp(argv[i],"-s")==0)
                strcpy(data_file, argv[i+1]);
            else if(strcmp(argv[i],"-o")==0)
                strcpy(out_file, argv[i+1]);
            else if(strcmp(argv[i],"-c")==0)
                strcpy(solution_file, argv[i+1]);
            else if(strcmp(argv[i],"-t")==0)
                MAX_TIME = atoi(argv[i+1]);
        }
        //printf("data_file= %s, output_file= %s, sln_file=%s, max_time=%d", data_file, out_file, solution_file, MAX_TIME);
    }
    struct problem_struct** my_problems = load_problems(data_file);

    if(strlen(solution_file)>=0)
    {
        if(strcmp(out_file,"")==0) strcpy(out_file, "my_solutions.txt"); //default output
        FILE* pfile = fopen(out_file, "w"); //open a new file
        fprintf(pfile, "%d\n", num_of_problems); fclose(pfile);
        for(int k=0; k<num_of_problems; k++)
        {
            best_sln.objective=0; best_sln.feasibility=0;
            for(int run=0; run<NUM_OF_RUNS; run++)
            {
                srand(RAND_SEED[run]);
                switch (alg)
                {
                case SA:
                    printf("Running Simulated Annealing...\n");
                    SimulatedAnnealing(my_problems[k]); // call SA method
                    break;
                case TS:
                    printf("Running Tabu Search...\n");
                    TabuSearch(my_problems[k]);
                    break;
                case VNS:
                    printf("Running VNS...\n");
                    varaible_neighbourhood_search(my_problems[k]);
                    break;
                default:
                    printf("No algorithm selected. Please select an algorithm!\n");
                    return 1;
                }
            }
            output_solution(&best_sln,out_file);
        }
    }
    for(int k=0; k<num_of_problems; k++)
    {
       free_problem(my_problems[k]); //free problem data memory
    }
    free(my_problems); //free problems array
    if(best_sln.x!=NULL && best_sln.cap_left!=NULL){ free(best_sln.cap_left); free(best_sln.x);} //free global
    return 0;
}
