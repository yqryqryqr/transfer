#include <iostream>
#include <vector>
#include <cmath>
#include <float.h>
#include <fstream>

struct pvar{
    double rho, v, p;
};

struct cvar{
    double rho, rhov, E;
};

double Gamma = 1.4; 
double CFL = 0.8;   
double x0 = 0, x1 = 1; 
double T_end = 0.25;  
int nCells = 100; 
int nPoints = nCells + 1; 
int ghost_cells = 3;
double dx = (x1 - x0)/(nPoints - 2);   

void initializeGrid(std::vector<pvar> & grid);
pvar con2pri(double rho, double rhov, double E, double Gamma);
cvar pri2con(double rho, double v, double p, double Gamma);
std::vector<double> f(double rho, double v, double p, double E);
double Minbee(double r);
std::vector<double> fri_flux(const cvar& uL, const cvar& uR, double pL, double pR, double dt, double dx);
std::vector<double> RI_flux(const cvar& uL, const cvar uR, double pL, double pR, double dt, double dx, double Gamma);
std::vector<double> FORCE(const cvar& uL, const cvar& uR, double pL, double pR, double dt, double dx, double Gamma);
void data_reconstruction(const std::vector<double> u, std::vector<double>& uL, std::vector<double>& uR);
void half_time(cvar &uL_half, cvar &uR_half, cvar &uL, cvar& uR, double dt, double dx, double Gamma);
double calcs(double p, double rho, double Gamma);
double update(const std::vector<cvar>& grid, double dx, double Gamma, double CFL);

void initializeGrid(std::vector<pvar>& grid) {
    for(int i = ghost_cells; i < nCells + ghost_cells; i++){
        double x = (i - ghost_cells + 0.5)* dx;
        if (x < 0.5){
            grid[i] = {1.0, 0, 1};
        }else{
            grid[i] = {0.125, 0.0, 0.1};
        }
    }
}

pvar con2pri(double rho, double rhov, double E, double Gamma){
    pvar pri;
    pri.rho = rho;
    pri.v = rhov / rho;
    pri.p = (Gamma - 1) * (E - 0.5 * rhov * rhov / rho);
    return pri;
}

cvar pri2con(double rho, double v, double p , double Gamma){
    cvar con;
    con.rho = rho;
    con.rhov = rho * v;
    con.E = p / (Gamma - 1) + 0.5 * rho * v * v; //E = p/(Gamma-1) + 0.5 rho v^2
    return con;
}

//调用方式
// 从原始变量转换为守恒变量
//pvar pri = {rho, v, p};  // 假设有原始变量
//cvar con = pri2con(pri.rho, pri.v, pri.p, Gamma);

// 访问计算后的守恒变量
//double rho_con = con.rho;
//double rhov_con = con.rhov;
//double E_con = con.E;

std::vector<double> f(double rho, double v, double p, double E){
    return{rho * v, rho * v * v + p, (E + p) * v};
}

double Minbee(double r){
    if(r <= 0) return 0;
    if(r <= 1) return r;
    return std::min(1.0, 2.0 / (1 + r));
}

std::vector<double> fri_flux(const cvar& uL, const cvar& uR, double pL, double pR, double dt, double dx){
    std::vector<double> flux_f(3);
    //计算左右通量
    std::vector<double> FL = f(uL.rho, uL.rhov / uL.rho, pL, uL.E);
    std::vector<double> FR = f(uR.rho, uR.rhov / uR.rho, pR, uR.E);
    flux_f[0] = 0.5 * (FL[0] + FR[0]) - 0.5 * (dx / dt) * (uR.rho - uL.rho);
    flux_f[1] = 0.5 * (FL[1] + FR[1]) - 0.5 * (dx / dt) * (uR.rhov - uL.rhov);
    flux_f[2] = 0.5 * (FL[2] + FR[2]) - 0.5 * (dx / dt) * (uR.E - uL.E);
    return flux_f;
}
//调用方式：
//cvar U_L = pri2con(rho_L, v_L, p_L, Gamma);
//cvar U_R = pri2con(rho_R, v_R, p_R, Gamma);
//std::vector<double> flux = laxFriedrichsFlux(U_L, U_R, p_L, p_R, dt, dx);

std::vector<double> RI_flux(const cvar& uL, const cvar uR, double pL, double pR, double dt, double dx, double Gamma){
    std::vector<double> flux(3);
    std::vector<double> FL = f(uL.rho, uL.rhov / uL.rho, pL, uL.E);
    std::vector<double> FR = f(uR.rho, uR.rhov / uR.rho, pR, uR.E);
    //计算半实践部的中间状态u*
    cvar u_star;
    u_star.rho = 0.5 * (uL.rho + uR.rho) - 0.5 * (dt/dx) * (FR[0] - FL[0]);
    u_star.rhov = 0.5 * (uL.rhov + uR.rhov) - 0.5 * (dt/dx) * (FR[1] - FL[1]);
    u_star.E = 0.5 * (uL.E + uR.E) - 0.5 * (dt/dx) * (FR[2] - FL[2]);

    //计算u*的通量
    pvar pri_star = con2pri(u_star.rho, u_star.rhov, u_star.E, Gamma);
    std::vector<double> flux_r = f(pri_star.rho, pri_star.v, pri_star.p, u_star.E);
    return flux_r;
}

std::vector<double> FORCE(const cvar& uL, const cvar& uR, double pL, double pR, double dt, double dx, double Gamma){
    std::vector<double> flux(3);
    std::vector<double> flux_f = fri_flux(uL, uR, pL, pR, dt, dx);
    std::vector<double> flux_r = RI_flux(uL, uR, pL, pR, dt, dx, Gamma);
    for(int i = 0; i < 3; i++){
        flux[i] = 0.5 * (flux_f[i] + flux_r[i]);
    }
    return flux;
}

void data_reconstruction(double uL, double u, double uR, double &uL_bar, double &uR_bar){
    double delta0, delta1, delta, r;
    delta0 = u - uL;//向后
    delta1 = uR- u;//向前
    delta = 0.5 * (delta0 + delta1);
    if (delta1 == 0){
        r = 0;
    }else{
        r = delta0 / delta1;
    }
    double phi;
    phi = Minbee(r);

    uL_bar = u - 0.5 * phi * delta;
    uR_bar = u + 0.5 * phi * delta;
}

void half_time(cvar &uL_half, cvar &uR_half, cvar &uL, cvar& uR, double dt, double dx, double Gamma){

    pvar priL = con2pri(uL.rho, uL.rhov, uL.E, Gamma);
    pvar priR = con2pri(uR.rho, uR.rhov, uR.E, Gamma);

    std::vector<double> FL = f(priL.rho, priL.v, priL.p, uL.E);
    std::vector<double> FR = f(priR.rho, priR.v, priR.p, uR.E);
    
    uL_half.rho = uL.rho - 0.5 * (dt/dx) * (FR[0] - FL[0]);
    uL_half.rhov = uL.rhov - 0.5 * (dt/dx) * (FR[1] - FL[1]);
    uL_half.E = uL.E - 0.5 * (dt/dx) * (FR[2] - FL[2]);

    uR_half.rho = uR.rho - 0.5 * (dt/dx) * (FR[0] - FL[0]);
    uR_half.rhov = uR.rhov- 0.5 * (dt/dx) * (FR[1] - FL[1]);
    uR_half.E = uR.E - 0.5 * (dt/dx) * (FR[2] - FL[2]);
}

//std::vector<double> SLIC(std::vector<cvar>& grid, double dt, double dx, double Gamma)
void SLIC(std::vector<std::vector<double>>& flux, const std::vector<cvar>& grid, double dt, double dx, double Gamma){
    
    std::vector<cvar> halfL(nCells); 
    std::vector<cvar> halfR(nCells);
    
    for(int i = ghost_cells - 1; i < nCells + ghost_cells + 1; i++){
        cvar uL = grid[i - 1];
        cvar u = grid[i];
        cvar uR = grid[i + 1];
        
        double rho, v, p, E;
        
        pvar priL = con2pri(uL.rho, uL.rhov, uL.E, Gamma);
        pvar pri = con2pri(u.rho, u.rhov, u.E, Gamma);
        pvar priR = con2pri(uR.rho, uR.rhov, uR.E, Gamma);
        double rhoL_bar, rhoR_bar, vL_bar, vR_bar, pL_bar, pR_bar;
    
        data_reconstruction(priL.rho, pri.rho, priR.rho, rhoL_bar,rhoR_bar);
        data_reconstruction(priL.v, pri.v, priR.v, vL_bar,vR_bar);
        data_reconstruction(priL.p, pri.p, priR.p, pL_bar,pR_bar);
         // 检查 NaN 传播
         if (std::isnan(pri.p) || std::isnan(pri.rho) || std::isnan(pri.v)) {
            std::cerr << "ERROR: NaN detected in data reconstruction at i=" << i << std::endl;
            exit(1);
        }

        // std::cout <<"data reconstruction done at i=" << i << std::endl;
    
        cvar uL_half, uR_half;
        cvar uL_bar, uR_bar;
    
        uL_bar = pri2con(rhoL_bar, vL_bar, pL_bar, Gamma);
        uR_bar = pri2con(rhoR_bar, vR_bar, pR_bar, Gamma);
        
        half_time(uL_half, uR_half, uL_bar, uR_bar, dt, dx, Gamma);
        // std::cout <<"half_time done" << std::endl;

    
        halfL[i] = uL_half;
        halfR[i] = uR_half;
    }
    
    for(int i = ghost_cells - 1; i < nCells + ghost_cells; i++){
        cvar uL = halfR[i];
        cvar uR = halfL[i + 1];
        pvar priL = con2pri(uL.rho, uL.rhov, uL.E, Gamma);
        pvar priR = con2pri(uR.rho, uR.rhov, uR.E, Gamma);
        // std::cout << "convert done" << std::endl;
        flux[i] = FORCE(uL, uR, priL.p, priR.p, dt, dx, Gamma);
        // std::cout << "FORCE done" << std::endl;

    }
}

double calcs(double p, double rho, double Gamma){
    return sqrt(Gamma * p / rho);
}

double update(const std::vector<cvar>& grid, double dx, double Gamma, double CFL){
    double max_speed = 0.0;
    for(int i = 0; i < grid.size(); i++){
        pvar pri = con2pri(grid[i].rho, grid[i].rhov, grid[i].E, Gamma);
        double cs = calcs(pri.p, pri.rho, Gamma);
        double speed = fabs(pri.v) + cs;
        if (speed > max_speed){
            max_speed = speed;
        }
    }
    double dt =  CFL * dx / max_speed;

    // 打印 dt 和 max_speed，检查是否逐渐趋近于 0
    std::cout << "Max speed: " << max_speed << ", Computed dt: " << dt << std::endl;

    // if (dt < 1e-10) {
    //     std::cerr << "ERROR: dt too small, possible numerical instability!" << std::endl;
    //     exit(1);
    // }
}

int main(){
    // std::cout << "simulation started." << std::endl;
    std::vector<pvar> primitive_grid(nCells + 2 * ghost_cells);
    // std::cout << "initialise start" << std::endl;
    initializeGrid(primitive_grid);
    for(int i = 0; i < ghost_cells; i++){
        primitive_grid[i] = primitive_grid[ghost_cells];
        primitive_grid[nCells + ghost_cells + i] = primitive_grid[nCells + ghost_cells - 1];
    }
    
    std::vector<cvar> grid(nCells + 2 * ghost_cells);
    std::cout <<"converting to conserved variables" << std::endl;
    for (int i = 0; i < nCells + 2 * ghost_cells; i++){
        grid[i] = pri2con(primitive_grid[i].rho, primitive_grid[i].v, primitive_grid[i].p, Gamma);
    }

    // 计算初始时间步长
    
    double t = 0.0;
    // do-while
    do {
        std::cout << "loop start" << std::endl;
        double dt = update(grid, dx, Gamma, CFL);
        std::cout << "dt calculated, entering SLIC()" << std::endl;

        std::vector<std::vector<double>> flux(nCells + 2 * ghost_cells, std::vector<double>(3));
        
        SLIC(flux, grid, dt, dx, Gamma);
        std::cout << "SLIC done" << std::endl;
        std::vector<cvar> new_grid = grid;
        for (int i = ghost_cells; i < nCells + ghost_cells; i++) {  
            new_grid[i].rho  -= (dt / dx) * (flux[i][0] - flux[i - 1][0]);
            new_grid[i].rhov -= (dt / dx) * (flux[i][1] - flux[i - 1][1]);
            new_grid[i].E    -= (dt / dx) * (flux[i][2] - flux[i - 1][2]);
        }
        
        //
        for(int i = 0; i < ghost_cells; i++){
            new_grid[i] = new_grid[ghost_cells];
            new_grid[nCells + ghost_cells + i] = new_grid[nCells + ghost_cells - 1];
        }

        // 更新网格
        grid = new_grid;

//      dt = update(grid, dx, Gamma, CFL);
        t += dt;

        // 输出进度
        
    }while (t < T_end);

    std::ofstream outFile("SLIC_results.dat");
    for (int i = ghost_cells; i < nCells + ghost_cells; i++) {
        pvar pri = con2pri(grid[i].rho, grid[i].rhov, grid[i].E, Gamma);
        outFile << i * dx << " " << pri.rho << " " << pri.v << " " << pri.p << std::endl;
    }
    outFile.close();

    std::cout << "Simulation completed using FORCE method. Results saved in 'force_results.txt'." << std::endl;
    return 0;
}
