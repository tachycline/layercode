#from matplotlib import pyplot as plt
#from matplotlib import animation

#    def plot_psi(self):
#        return plt.contour(self.xx, self.yy, self.psi)
        
#    def plot_q(self):
#        return plt.contour(self.xx, self.yy, self.q)
        


#    def animate_results(self, filename, stride=1):
#        fig = plt.figure(figsize=(10,10))
#        ax = plt.axes()
#
#        plt.xlabel(r'x')
#        plt.ylabel(r'y')
#
#        def animate(i):
#            qhat = self.results[int(i*stride)]
#            psihat = self.get_psihat_from_qhat(qhat)
#            z = np.fft.irfft2(psihat)
#            ax.clear()
#            cont = plt.contour(self.xx, self.yy, z)
#            return cont
#
#        anim = animation.FuncAnimation(fig, animate, frames=len(self.results)//stride, blit=False)
#        mywriter = animation.FFMpegWriter()
#        anim.save(filename, bitrate=10000)
