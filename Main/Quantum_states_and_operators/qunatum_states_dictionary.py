"Bulid all the density matrix states"

class rhoMatrixNames:

    def __init__(self,N):
        self.N = N
        self.rho = {}
        counter = 0
        for i in range(N):
            for j in range(N):
                state = 'rho' + str(i + 1) + str(j + 1)
                self.rho.update({state: counter})
                counter += 1

    def getLocationByName(self,rho_name):
        return self.rho[rho_name]

    @staticmethod
    def getStatesNames(N):
        rho={}
        counter=0
        for i in range(N):
            for j in range(N):
                state='rho'+str(i+1)+str(j+1)
                rho.update({state:counter})
                counter+=1
        return rho

    @staticmethod
    def buildRhoName(row,col):
        ret_val = 'rho'+str(row)+str(col)
        return ret_val


    def __str__(self):
        txt = ''
        for i in range(self.N):
            for idx, j in enumerate(range(self.N-1)):
                txt +='rho'+str(i+1)+str(j+1)+'| '
            txt += 'rho' + str(i + 1) + str(j + 2)
            txt += '\n'
            txt += 30 * '-'
            txt += '\n' 
        return txt

