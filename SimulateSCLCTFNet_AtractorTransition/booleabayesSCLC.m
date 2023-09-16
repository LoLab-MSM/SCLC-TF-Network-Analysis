function [networkStates, trackReaction, transitionStatus, message] = booleabayesSCLC(numIter, initialState, rules, attractors)

    %% Simulation Settings:
    networkStates = zeros(35,numIter+1);
    networkStates(:,1) = initialState;
    trackReaction = [];
    transitionStatus = 0;
    message = sprintf('Transition have not been observed within %d iterations', numIter);

    %% Assign the lookup tables:
    ISL1_table = rules{1}; REST_table = rules{2}; STAT6_table = rules{3};
    TEAD4_table = rules{4}; ZNF217_table = rules{5}; NEUROD2_table = rules{6}; 
    HES1_table = rules{7}; SIX5_table = rules{8}; BCL3_table = rules{9};
    CEBPD_table = rules{10}; SMAD4_table = rules{11}; FOXA2_table = rules{12};
    ELF3_table = rules{13}; NR0B1_table = rules{14}; RARG_table = rules{15};
    MYCN_table = rules{16}; RBPJ_table = rules{17}; NR0B2_table = rules{18};
    FOXA1_table = rules{19}; FLI1_table = rules{20}; NEUROD1_table = rules{21};
    GATA4_table = rules{22}; OLIG2_table = rules{23}; ASCL1_table = rules{24}; 
    GFI1B_table = rules{25}; RCOR2_table = rules{26}; SOX11_table = rules{27};
    POU2F3_table = rules{28}; MITF_table = rules{29}; TCF3_table = rules{30};
    TCF4_table = rules{31}; YAP1_table = rules{32}; KLF2_table = rules{33};
    MYC_table = rules{34}; EBF1_table = rules{35};

    %% RunSimulation
    for i = 2:numIter+1

        ASCL1 = networkStates(1,i-1); FOXA1 = networkStates(2,i-1); FOXA2 = networkStates(3,i-1);
        ELF3 = networkStates(4,i-1); RBPJ = networkStates(5,i-1); FLI1 = networkStates(6,i-1);
        SMAD4 = networkStates(7,i-1); NR0B2 = networkStates(8,i-1); NR0B1 = networkStates(9,i-1);
        BCL3 = networkStates(10,i-1); STAT6 = networkStates(11,i-1); ISL1 = networkStates(12,i-1);
        SOX11 = networkStates(13,i-1); CEBPD = networkStates(14,i-1); EBF1 = networkStates(15,i-1);
        TCF4 = networkStates(16,i-1); RCOR2 = networkStates(17,i-1); TCF3 = networkStates(18,i-1);
        NEUROD2 = networkStates(19,i-1); OLIG2 = networkStates(20,i-1); MITF = networkStates(21,i-1);
        SIX5 = networkStates(22,i-1); TEAD4 = networkStates(23,i-1); ZNF217 = networkStates(24,i-1);
        KLF2 = networkStates(25,i-1); GATA4 = networkStates(26,i-1); REST = networkStates(27,i-1);
        HES1 = networkStates(28,i-1); RARG = networkStates(29,i-1); MYCN = networkStates(30,i-1);
        NEUROD1 = networkStates(31,i-1); GFI1B = networkStates(32,i-1); POU2F3 = networkStates(33,i-1);
        YAP1 = networkStates(34,i-1); MYC = networkStates(35,i-1);

        reactID = randi([1, 35], 1);
        trackReaction = [trackReaction, reactID];


        % ASCL1:
%         if(reactID == 1)
%             probASCL1On = ASCL1_table(1, getIndex([FLI1 TEAD4 OLIG2 HES1 MYC MITF SMAD4 KLF2]));
%             if(rand() <= probASCL1On)
%                 networkStates(1,i) = 1;
%             end
%         else
%             networkStates(1,i) = networkStates(1,i-1);
%         end
        networkStates(1,i) = 1; % ON for FLI1-ASCL1-MITF pathway

        % FOXA1:
        if(reactID == 2)
            probFOXA1On = FOXA1_table(1, getIndex([ASCL1 NEUROD2 ELF3 NR0B1 FOXA2 FLI1 YAP1 ISL1 MYCN TCF3 OLIG2 ZNF217 MYC]));
            if(rand() <= probFOXA1On)
                networkStates(2,i) = 1;
            end
        else
            networkStates(2,i) = networkStates(2,i-1);
        end

        % FOXA2:
        if(reactID == 3)
            probFOXA2On = FOXA2_table(1, getIndex([ASCL1 ISL1 FLI1 REST TEAD4 TCF3 FOXA1 SMAD4 EBF1 OLIG2]));
            if(rand() <= probFOXA2On)
                networkStates(3,i) = 1;
            end
        else
            networkStates(3,i) = networkStates(3,i-1);
        end

        % ELF3:
        if(reactID == 4)
            probELF3On = ELF3_table(1, getIndex([FOXA1 GATA4 MYC ZNF217 MITF EBF1 RBPJ]));
            if(rand() <= probELF3On)
                networkStates(4,i) = 1;
            end
        else
            networkStates(4,i) = networkStates(4,i-1);
        end

        % RBPJ:
        if(reactID == 5)
            probRBPJOn = RBPJ_table(1, getIndex([ASCL1 ISL1 FOXA2 GATA4 TCF4 EBF1 TCF3 OLIG2 MITF GFI1B SMAD4 RCOR2]));
            if(rand() <= probRBPJOn)
                networkStates(5,i) = 1;
            end
        else
            networkStates(5,i) = networkStates(5,i-1);
        end

        % FLI1:
%         if(reactID == 6)
%             probFLI1On = FLI1_table(1, getIndex([FOXA2 MYC SOX11 SMAD4 STAT6 FOXA1 TCF3]));
%             if(rand() <= probFLI1On)
%                 networkStates(6,i) = 1;
%             end
%         else
%             networkStates(6,i) = networkStates(6,i-1);
%         end
        networkStates(6,i) = 1; % ON for FLI1-ASCL1-MITF pathway

        % SMAD4:
        if(reactID == 7)
            probSMAD4On = SMAD4_table(1, getIndex([FOXA1 FOXA2 MYC NEUROD2]));
            if(rand() <= probSMAD4On)
                networkStates(7,i) = 1;
            end
        else
            networkStates(7,i) = networkStates(7,i-1);
        end

        % NR0B2:
        if(reactID == 8)
            probNR0B2On = NR0B2_table(1, getIndex([ASCL1 FOXA2 FLI1 FOXA1 MYC NEUROD2 GATA4 RBPJ HES1]));
            if(rand() <= probNR0B2On)
                networkStates(8,i) = 1;
            end
        else
            networkStates(8,i) = networkStates(8,i-1);
        end

        % NR0B1:
        if(reactID == 9)
            probNR0B1On = NR0B1_table(1, getIndex(NR0B1)); % !!!
            if(rand() <= probNR0B1On)
                networkStates(9,i) = 1;
            end
        else
            networkStates(9,i) = networkStates(9,i-1);
        end

        % BCL3:
        if(reactID == 10)
            probBCL3On = BCL3_table(1, getIndex([SOX11 ELF3 REST MYC NR0B1]));
            if(rand() <= probBCL3On)
                networkStates(10,i) = 1;
            end
        else
            networkStates(10,i) = networkStates(10,i-1);
        end

        % STAT6:
        if(reactID == 11)
            probSTAT6On = STAT6_table(1, getIndex([TCF4 BCL3 FOXA2 GATA4 FLI1 MYC FOXA1 CEBPD MITF]));
            if(rand() <= probSTAT6On)
                networkStates(11,i) = 1;
            end
        else
            networkStates(11,i) = networkStates(11,i-1);
        end

        % ISL1:
        if(reactID == 12)
            probISL1On = ISL1_table(1, getIndex([FLI1 SMAD4 MYC SOX11 TCF4 GATA4 NR0B1]));
            if(rand() <= probISL1On)
                networkStates(12,i) = 1;
            end
        else
            networkStates(12,i) = networkStates(12,i-1);
        end

        % SOX11:
        if(reactID == 13)
            probSOX11On = SOX11_table(1, getIndex([ISL1 FLI1 TCF4 RBPJ MYC TCF3 FOXA1 OLIG2]));
            if(rand() <= probSOX11On)
                networkStates(13,i) = 1;
            end
        else
            networkStates(13,i) = networkStates(13,i-1);
        end

        % CEBPD:
        if(reactID == 14)
            probCEBPDOn = CEBPD_table(1, getIndex(CEBPD)); % !!!
            if(rand() <= probCEBPDOn)
                networkStates(14,i) = 1;
            end
        else
            networkStates(14,i) = networkStates(14,i-1);
        end

        % EBF1:
        if(reactID == 15)
            probEBF1On = EBF1_table(1, getIndex([YAP1 SMAD4 REST TEAD4 MYC CEBPD FLI1 TCF4 MITF]));
            if(rand() <= probEBF1On)
                networkStates(15,i) = 1;
            end
        else
            networkStates(15,i) = networkStates(15,i-1);
        end

        % TCF4:
        if(reactID == 16)
            probTCF4On = TCF4_table(1, getIndex([NEUROD2 ELF3 REST SMAD4 TEAD4 MYC GFI1B FOXA1 TCF3 MITF]));
            if(rand() <= probTCF4On)
                networkStates(16,i) = 1;
            end
        else
            networkStates(16,i) = networkStates(16,i-1);
        end

        % RCOR2:
        if(reactID == 17)
            probRCOR2On = RCOR2_table(1, getIndex([NEUROD2 OLIG2 FLI1 NR0B1 GATA4 TCF3 MYCN GFI1B TCF4 MYC]));
            if(rand() <= probRCOR2On)
                networkStates(17,i) = 1;
            end
        else
            networkStates(17,i) = networkStates(17,i-1);
        end

        % TCF3:
        if(reactID == 18)
            probTCF3On = TCF3_table(1, getIndex([STAT6 ISL1 MYC NR0B1 FLI1 YAP1 TCF4 MYCN]));
            if(rand() <= probTCF3On)
                networkStates(18,i) = 1;
            end
        else
            networkStates(18,i) = networkStates(18,i-1);
        end

        % NEUROD2:
        if(reactID == 19)
            probNEUROD2On = NEUROD2_table(1, getIndex([RCOR2 FLI1 NR0B2 REST MYC]));
            if(rand() <= probNEUROD2On)
                networkStates(19,i) = 1;
            end
        else
            networkStates(19,i) = networkStates(19,i-1);
        end

        % OLIG2:
        if(reactID == 20)
            probOLIG2On = OLIG2_table(1, getIndex([ASCL1 FLI1 NEUROD2 NR0B1 FOXA1 SMAD4 TCF4 TCF3]));
            if(rand() <= probOLIG2On)
                networkStates(20,i) = 1;
            end
        else
            networkStates(20,i) = networkStates(20,i-1);
        end

        % MITF:
        if(reactID == 21)
            probMITFOn = MITF_table(1, getIndex([SMAD4 YAP1 MYC SIX5 OLIG2 ZNF217 TCF4]));
            if(rand() <= probMITFOn)
                networkStates(21,i) = 1;
            end
        else
            networkStates(21,i) = networkStates(21,i-1);
        end

        % SIX5:
        if(reactID == 22)
            probSIX5On = SIX5_table(1, getIndex([MYC ELF3 BCL3 FOXA1 MITF]));
            if(rand() <= probSIX5On)
                networkStates(22,i) = 1;
            end
        else
            networkStates(22,i) = networkStates(22,i-1);
        end

        % TEAD4:
        if(reactID == 23)
            probTEAD4On = TEAD4_table(1, getIndex([ELF3 MYC CEBPD SMAD4 YAP1 FOXA1 MITF GATA4 MYCN]));
            if(rand() <= probTEAD4On)
                networkStates(23,i) = 1;
            end
        else
            networkStates(23,i) = networkStates(23,i-1);
        end

        % ZNF217:
        if(reactID == 24)
            probZNF217On = ZNF217_table(1, getIndex(FOXA2));
            if(rand() <= probZNF217On)
                networkStates(24,i) = 1;
            end
        else
            networkStates(24,i) = networkStates(24,i-1);
        end

        % KLF2:
        if(reactID == 25)
            probKLF2On = KLF2_table(1, getIndex([SMAD4 NR0B1 CEBPD FLI1 YAP1 TCF3 TCF4 MYC]));
            if(rand() <= probKLF2On)
                networkStates(25,i) = 1;
            end
        else
            networkStates(25,i) = networkStates(25,i-1);
        end

        % GATA4:
        if(reactID == 26)
            probGATA4On = GATA4_table(1, getIndex([FLI1 SMAD4 MYC TCF4 FOXA2 TEAD4 YAP1]));
            if(rand() <= probGATA4On)
                networkStates(26,i) = 1;
            end
        else
            networkStates(26,i) = networkStates(26,i-1);
        end

        % REST:
        if(reactID == 27)
            probRESTOn = REST_table(1, getIndex([YAP1 NR0B1 MYC GATA4 TCF4 FOXA1 TCF3]));
            if(rand() <= probRESTOn)
                networkStates(27,i) = 1;
            end
        else
            networkStates(27,i) = networkStates(27,i-1);
        end

        % HES1:
        if(reactID == 28)
            probHES1On = HES1_table(1, getIndex([NEUROD2 GATA4 ASCL1 NR0B1 ISL1 FLI1 GFI1B YAP1 STAT6 OLIG2 MYCN TCF3 SMAD4 REST FOXA2 FOXA1 MYC ZNF217 MITF]));
            if(rand() <= probHES1On)
                networkStates(28,i) = 1;
            end
        else
            networkStates(28,i) = networkStates(28,i-1);
        end

        % RARG:
        if(reactID == 29)
            probRARGOn = RARG_table(1, getIndex([FLI1 NR0B1 SOX11 MYC FOXA1 NEUROD2 MYCN TCF3]));
            if(rand() <= probRARGOn)
                networkStates(29,i) = 1;
            end
        else
            networkStates(29,i) = networkStates(29,i-1);
        end

        % MYCN:
        if(reactID == 30)
            probMYCNOn = MYCN_table(1, getIndex([NEUROD2 YAP1 OLIG2 NR0B1 REST FLI1 RCOR2 TEAD4 TCF3 TCF4 FOXA1 RBPJ MYC]));
            if(rand() <= probMYCNOn)
                networkStates(30,i) = 1;
            end
        else
            networkStates(30,i) = networkStates(30,i-1);
        end

        % NEUROD1:
        if(reactID == 31)
            probNEUROD1On = NEUROD1_table(1, getIndex([SOX11 FLI1 EBF1 NR0B2 NR0B1 REST TCF3 HES1 ZNF217 OLIG2 SMAD4]));
            if(rand() <= probNEUROD1On)
                networkStates(31,i) = 1;
            end
        else
            networkStates(31,i) = networkStates(31,i-1);
        end

        % GFI1B:
        if(reactID == 32)
            probGFI1BOn = GFI1B_table(1, getIndex([POU2F3 TCF4 EBF1 REST FLI1]));
            if(rand() <= probGFI1BOn)
                networkStates(32,i) = 1;
            end
        else
            networkStates(32,i) = networkStates(32,i-1);
        end

        % POU2F3:
        if(reactID == 33)
            probPOU2F3On = POU2F3_table(1, getIndex(POU2F3));
            if(rand() <= probPOU2F3On)
                networkStates(33,i) = 1;
            end
        else
            networkStates(33,i) = networkStates(33,i-1);
        end

        % YAP1:
        if(reactID == 34)
            probYAP1On = YAP1_table(1, getIndex([FOXA2 REST ELF3 GATA4 FLI1 MITF SMAD4 FOXA1 RBPJ RARG MYC]));
            if(rand() <= probYAP1On)
                networkStates(34,i) = 1;
            end
        else
            networkStates(34,i) = networkStates(34,i-1);
        end

        % MYC:
        if(reactID == 35)
            probMYCOn = MYC_table(1, getIndex([FLI1 GATA4 ISL1 OLIG2 YAP1 NEUROD1 TCF4 SMAD4 CEBPD TCF3 FOXA1 ZNF217]));
            if(rand() <= probMYCOn)
                networkStates(35,i) = 1;
            end
        else
            networkStates(35,i) = networkStates(35,i-1);
        end
        
        % Check if a transition occur from Non-NE to another subtype"
        if((networkStates(1:27,i) == attractors(1,:)') | (networkStates(1:27,i) == attractors(2,:)'))
            message = 'Transition occured from Non-NE to NE subtype!';
            transitionStatus = 1;
            break
        elseif((networkStates(1:27,i) == attractors(3,:)') | (networkStates(1:27,i) == attractors(4,:)'))
            message = 'Transition occured from Non-NE to NEv1 subtype!';
            transitionStatus = 1;
            break
        elseif((networkStates(1:27,i) == attractors(5,:)') | (networkStates(1:27,i) == attractors(6,:)') | (networkStates(1:27,i) == attractors(7,:)') | (networkStates(1:27,i) == attractors(8,:)'))
            message = 'Transition occured from Non-NE to NEv2 subtype!';
            transitionStatus = 1;
            break
        end


    end

    %% Convert upstream binary combinations into a single index to retrieve its
    % associated probability from its' rule table.
    function index = getIndex(upstream)

        index = bin2dec(num2str(upstream)) + 1;

    end

end