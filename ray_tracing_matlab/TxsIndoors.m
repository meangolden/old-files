function [Txs] = TxsIndoors(X,Y,Z,fc,txAntenna,txPower)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
Txs = txsite() ; i=1;
for x = X
    for y = Y
        for z = Z
            t = txsite("cartesian", ...
            "Antenna",txAntenna, ...
            "AntennaPosition",[x;y;z], ...
            'TransmitterFrequency',fc, ...
            "TransmitterPower", txPower);
            Txs(i)=t;
            i=i+1;
        end
    end
end
