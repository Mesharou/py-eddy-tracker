from datetime import datetime, timedelta
from py_eddy_tracker import start_logger
from py_eddy_tracker.dataset.grid import UnRegularGridDataset
from numpy import zeros, arange, float
from netCDF4 import Dataset

# for plotting
from matplotlib import pyplot as plt

import sys,os,shutil


class RomsDataset(UnRegularGridDataset):
    
    __slots__ = list()
    
    def init_speed_coef(self, uname, vname):
        # xi_u and eta_v must be specified because this dimension are not use in lon/lat
        print(self.indexs)
        u = self.grid(uname, indexs=self.indexs)
        v = self.grid(vname, indexs=self.indexs)
        print('u.shape',u.shape)

        u = self.rho_2d(u.T).T
        v = self.rho_2d(v)
        self._speed_norm = (v ** 2 + u ** 2) ** .5


    @staticmethod
    def rho_2d(x):
        """Transformation to have u or v on same grid than h
        """
        M, Lp = x.shape
        new_x = zeros((M + 1, Lp))
        new_x[1:-1] = .5 * (x[:-1] + x[1:])
        new_x[0] = new_x[1]
        new_x[-1] = new_x[-2]
        return new_x

    @staticmethod
    def psi2rho(var_psi):
        """Transformation to have vrt on same grid than h
        """
        M, L = var_psi.shape
        Mp=M+1; Lp=L+1
        Mm=M-1; Lm=L-1
        var_rho = zeros((Mp, Lp))
        var_rho[1:M,1:L]=0.25*(var_psi[0:Mm,0:Lm]+var_psi[0:Mm,1:L]+var_psi[1:M,0:Lm]+var_psi[1:M,1:L])
        var_rho[0,:]=var_rho[1,:]
        var_rho[Mp-1,:]=var_rho[M-1,:]
        var_rho[:,0]=var_rho[:,1]
        var_rho[:,Lp-1]=var_rho[:,L-1]
        return var_rho

##########################

varname = sys.argv[1]

print("varname id %s", varname)

#varname = 'zeta'

simul = 'gigatl6'

##########################

if __name__ == '__main__':
    start_logger().setLevel('DEBUG')



    ####
    # create case-specific folder
    folder = varname + '/' + sys.argv[2] + '/'

    try:
        os.mkdir(folder)
    except OSError:
        print ("Directory %s already exists" % folder)


    # copy script in folder
    shutil.copy('detection_gigatl6.py',folder)

    
    
    ####
    # define simulation
    
    if 'gigatl1' in simul:
        folder_simulation = '~/megatl/GIGATL1/GIGATL1_1h/SURF/final/'
    elif 'gigatl6' in simul:
        folder_simulation = './HIS/'
    
    # Pick a depth/isopycnal

    if varname == 'zeta':
        s_rho = -1
        depth = 0
    elif  varname == 'ow':
        #isopycnal
        filename = './iso/gigatl6_1h_isopycnal_section.01440.nc'
        nc = Dataset(filename,'r')
        isopycnals = nc.variables['isopycnal'][:]
        nc.close()

        s_rho =  1 #
        depth = isopycnals[s_rho]
        

    
    # Time loop
    #for time in range(34800, 34801):
    for time in range(1440, 6000):

        
        # hourly data
        #dtfile = 1
        # 12-hourly
        dtfile = 12

        tfile = 5*24//dtfile
        
        ###########
        realyear_origin = datetime(2004,1,15)
        date = realyear_origin + timedelta(days=float(time)*dtfile/24.)
            
        infiletime = time%tfile
        filetime = time - time%tfile
        
        
        # Using times:
        #filename = './GIGATL6/gigatl6_1h_horizontal_section.' + '{0:05}'.format(filetime) + '.nc'
        
        if varname =='ow':
            #filename = './GIGATL6/gigatl6_1h_horizontal_section.' + '{0:05}'.format(filetime) + '.nc'
            filename = './iso/gigatl6_1h_isopycnal_section.' + '{0:05}'.format(filetime) + '.nc'


        elif varname == 'zeta':

            
            # or using dates
            date1 = realyear_origin + timedelta(days=float(filetime)*dtfile/24.)
            date2 = date1 + timedelta(days=float(4.5))

            if 'gigatl1' in simul:
                filedate =  '{0:04}'.format(date1.year)+'-'+\
                            '{0:02}'.format(date1.month) + '-'+\
                            '{0:02}'.format(date1.day)

                filename = folder_simulation + 'gigatl1_surf.' + filedate + '.nc'
                gridname = '~/megatl/GIGATL1/gigatl1_grd.nc'
                
            elif 'gigatl6' in simul:
                filedate =  '{0:04}'.format(date1.year)+'-'+\
                            '{0:02}'.format(date1.month) + '-'+\
                            '{0:02}'.format(date1.day)+ '-'+\
                            '{0:04}'.format(date2.year)+'-'+\
                            '{0:02}'.format(date2.month) + '-' +\
                            '{0:02}'.format(date2.day)
                            
                print(filedate)
            
                filename = './HIS/GIGATL6_1h_inst_surf_' + filedate + '.nc'
                gridname = filename

        # Identification
        if varname=='zeta' and 'gigatl6' in simul:
            lon_name, lat_name = 'nav_lon_rho', 'nav_lat_rho'
        elif varname=='ow':
            lon_name, lat_name = 'lon', 'lat'
        else:
            lon_name, lat_name = 'lon_rho', 'lat_rho'
            
        # domain grid points: x1, x2, y1, y2
        x1, x2 = 400, 1200
        y1, y2 = 400, 1200

        x1, x2 = 0, 1600
        y1, y2 = 0, 2000

        h = RomsDataset(filename, lon_name, lat_name, gridname = gridname,
                indexs=dict(time=infiletime,time_counter=infiletime,
                eta_rho=slice(y1, y2),
                xi_rho=slice(x1, x2),
                eta_v=slice(y1, y2-1),
                xi_u=slice(x1, x2-1),
                y_rho=slice(y1, y2),
                x_rho=slice(x1, x2),
                y_v=slice(y1, y2-1),
                y_u=slice(y1, y2),
                x_u=slice(x1, x2-1),
                x_v=slice(x1, x2),
                s_rho=s_rho)
        )
        

        
        # Identification
        if varname=='zeta':
            z_min = -2 ; z_max = 1.5; step = 0.02
        elif varname=='ow':
            # ow is multiplied by 1e10
            z_min = -1; z_max = -0.01; step = 0.01


        a, c = h.eddy_identification(varname, 'u', 'v', date, z_min =  z_min, z_max = z_max, step = step, pixel_limit=(10, 2000), shape_error=40, force_height_unit='m',force_speed_unit='m/s',vorticity_name='vrt')
        
 
        ####

        filename = 'gigatl6_1h_' + varname + '_'+  '{0:04}'.format(depth) + '_'
        
        with Dataset(date.strftime(folder + 'Anticyclonic_' + filename + '%Y%m%d%H.nc'), 'w') as h:
            a.to_netcdf(h)
        with Dataset(date.strftime(folder + 'Cyclonic_' + filename + '%Y%m%d%H.nc'), 'w') as h:
            c.to_netcdf(h)
        
        # PLOT
        

        fig = plt.figure(figsize=(15,7))
        ax = fig.add_axes([.03,.03,.94,.94])
        ax.set_title('Eddies detected -- Cyclonic(red) and Anticyclonic(blue)')
        #ax.set_ylim(-75,75)
        ax.set_xlim(250,360)
        ax.set_aspect('equal')
        a.display(ax, color='b', linewidth=.5)
        c.display(ax, color='r', linewidth=.5)
        ax.grid()
        fig.savefig(folder + 'eddies_' + date.strftime( filename + '%Y%m%d%H') +'.png')


