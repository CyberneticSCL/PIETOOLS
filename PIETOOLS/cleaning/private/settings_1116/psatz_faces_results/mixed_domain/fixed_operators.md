# Fixed operators on the mixed spaces: positivity under the slack reader (FIX, next3.log, 11/16 02:30-03:05)

The constraint X >= 0 by lpi_ineq_sop for a fixed 'copvar' X, the certificate sized by the reader from the support of X (light tier: dD 1, dw 1, terms 'auto'), or with no terms, or with dD 2, dw 2; MOSEK. 'optimal' = X is in the slack's cone at that sizing. The batch was cut after S18_Pstar (every target certified; the raised-degree rows of P* and the derivative would only be larger).

| X | options | status | ndv | t (s) | D | w | Gram sizes N | eta |
|---|---|---|---|---|---|---|---|---|
| S17_TtT | {} | optimal | 14500 | 1 | [1 | [2 | [90 | 2.1e-08 |
| S17_TtT | {"psatz":"none"} | optimal | 8100 | 0 | [1 | [2 | 90 | 3.4e-09 |
| S17_TtT | {"dD":2,"dw":2} | optimal | 28420 | 0 | [2 | [2 | [126 | 1.8e-08 |
| S17_TtE1T | {} | optimal | 3625 | 0 | [1 | [2 | [45 | 3.1e-08 |
| S17_TtE1T | {"psatz":"none"} | optimal | 2025 | 0 | [1 | [2 | 45 | 2.1e-08 |
| S17_TtE1T | {"dD":2,"dw":2} | optimal | 7105 | 0 | [2 | [2 | [63 | 5.6e-09 |
| S17_TtE2T | {} | optimal | 14500 | 0 | [1 | [2 | [90 | 1.3e-08 |
| S17_TtE2T | {"psatz":"none"} | optimal | 8100 | 0 | [1 | [2 | 90 | 5.4e-09 |
| S17_TtE2T | {"dD":2,"dw":2} | optimal | 28420 | 0 | [2 | [2 | [126 | 4.5e-09 |
| S17_AtA | {} | optimal | 20880 | 0 | [1 | [2 | [108 | 1.3e-09 |
| S17_AtA | {"psatz":"none"} | optimal | 11664 | 0 | [1 | [2 | 108 | 1.1e-08 |
| S17_AtA | {"dD":2,"dw":2} | optimal | 37120 | 1 | [2 | [2 | [144 | 2.3e-08 |
| S3__TtT | {} | optimal | 255780 | 12 | [2 | [2 | [378 | 1.0e-08 |
| S3__TtT | {"psatz":"none"} | optimal | 142884 | 9 | [2 | [2 | 378 | 3.5e-09 |
| S3__TtT | {"dD":2,"dw":2} | optimal | 751680 | 61 | [3 | [2 | [648 | 1.1e-08 |
| S3__TtE1T | {} | optimal | 7105 | 0 | [2 | [2 | [63 | 1.0e-08 |
| S3__TtE1T | {"psatz":"none"} | optimal | 3969 | 0 | [2 | [2 | 63 | 3.8e-09 |
| S3__TtE1T | {"dD":2,"dw":2} | optimal | 11745 | 0 | [3 | [2 | [81 | 5.7e-09 |
| S3__TtE2T | {} | optimal | 255780 | 13 | [2 | [2 | [378 | 2.8e-09 |
| S3__TtE2T | {"psatz":"none"} | optimal | 142884 | 9 | [2 | [2 | 378 | 3.3e-09 |
| S3__TtE2T | {"dD":2,"dw":2} | optimal | 751680 | 69 | [3 | [2 | [648 | 1.2e-09 |
| S3__AtA | {} | optimal | 438625 | 18 | [2 | [2 | [495 | 2.7e-08 |
| S3__AtA | {"psatz":"none"} | optimal | 245025 | 14 | [2 | [2 | 495 | 3.7e-09 |
| S3__AtA | {"dD":2,"dw":2} | optimal | 1148545 | 104 | [3 | [2 | [801 | 1.1e-08 |
| S18_TtT | {} | optimal | 130500 | 4 | [2 | [2 | [270 | 1.3e-08 |
| S18_TtT | {"psatz":"none"} | optimal | 72900 | 3 | [2 | [2 | 270 | 7.0e-09 |
| S18_TtT | {"dD":2,"dw":2} | optimal | 454720 | 26 | [3 | [2 | [504 | 1.2e-08 |
| S18_TtE1T | {} | optimal | 7105 | 0 | [2 | [2 | [63 | 3.7e-09 |
| S18_TtE1T | {"psatz":"none"} | optimal | 3969 | 0 | [2 | [2 | 63 | 5.0e-09 |
| S18_TtE1T | {"dD":2,"dw":2} | optimal | 11745 | 0 | [3 | [2 | [81 | 6.7e-09 |
| S18_TtE2T | {} | optimal | 130500 | 4 | [2 | [2 | [270 | 1.4e-08 |
| S18_TtE2T | {"psatz":"none"} | optimal | 72900 | 3 | [2 | [2 | 270 | 1.6e-08 |
| S18_TtE2T | {"dD":2,"dw":2} | optimal | 454720 | 25 | [3 | [2 | [504 | 7.0e-09 |
| S18_AtA | {} | optimal | 198505 | 5 | [2 | [2 | [333 | 1.6e-08 |
| S18_AtA | {"psatz":"none"} | optimal | 110889 | 4 | [2 | [2 | 333 | 1.2e-08 |
| S18_AtA | {"dD":2,"dw":2} | optimal | 612625 | 31 | [3 | [2 | [585 | 3.7e-08 |
| S18_Pstar | {} | optimal | 4698000 | 776 | [3 | [3 | [1440 | 5.2e-07 |
