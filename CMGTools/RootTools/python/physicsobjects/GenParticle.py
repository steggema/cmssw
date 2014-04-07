from CMGTools.RootTools.physicsobjects.PhysicsObject import *

class GenParticle( PhysicsObject ):
    def __str__(self):
        base = super(GenParticle, self).__str__()
        theStr = '{base}, status = {status:>2}'.format(base=base, status=self.status())
        return theStr
    def getFinal(self):
        if self.numberOfDaughters() == 1 and self.daughter(0).pdgId() == self.pdgId():
            return GenParticle(self.daughter(0)).getFinal()
        return self



class GenLepton( GenParticle ):
    def sip3D(self):
        '''Just to make generic code work on GenParticles'''
        return 0
    def relIso(self, dummy):
        '''Just to make generic code work on GenParticles'''
        return 0

    def absIso(self, dummy):
        '''Just to make generic code work on GenParticles'''
        return 0

    def absEffAreaIso(self,rho):
        '''Just to make generic code work on GenParticles'''
        return 0

    def relEffAreaIso(self,rho):
        '''Just to make generic code work on GenParticles'''
        return 0
