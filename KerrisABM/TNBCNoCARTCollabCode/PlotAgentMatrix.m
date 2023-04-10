%PlotAgentMatrix

figure
            hold on
             p=patch(isosurface(AgentMatrix==1,0));
             set(p,'facecolor','red' ,'edgecolor', 'none');
                daspect([1 1 1])
                isonormals(AgentMatrix==1,p)
                view(3);

                camlight

                lighting gouraud 