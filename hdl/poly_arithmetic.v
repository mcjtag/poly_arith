`timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
// Company: 
// Engineer: Dmitry Matyunin (https://github.com/mcjtag)
// 
// Create Date: 10.08.2023 19:15:02
// Design Name: 
// Module Name: 
// Project Name: poly_arithmetic
// Target Devices:
// Tool Versions:
// Description: Polynomial arithmetic implementation
// Dependencies:
// Revision:
// Revision 0.01 - File Created
// Additional Comments:
// License: MIT
//  Copyright (c) 2023 Dmitry Matyunin
//  Permission is hereby granted, free of charge, to any person obtaining a copy
//  of this software and associated documentation files (the "Software"), to deal
//  in the Software without restriction, including without limitation the rights
//  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
//  copies of the Software, and to permit persons to whom the Software is
//  furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included in
//  all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
//  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
//  THE SOFTWARE.
// 
//////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////
//
//  Polynomial representation:
//
//  1) GF(2)
//	 if: 
//      p(x) = a[0] + a[1]*x + a[2]*x^2 + ... + a[k]*x^k (a[i] is in GF(2), k <= POLY_P_WIDTH-1)
//   then:
//      px[POLY_X_WIDTH-1:0] = {a[k], ..., a[2], a[1], a[0]}
//   
//  2) GF(2^m)
//   if:
//      p(x) = b[0] + b[1]*x + b[2]*x^2 + ... + b[k]*x^k (b[i] is in GF(2^m), k <= POLY_P_WIDTH-1)
//		b[i] = alpha^s, were s is from set: -inf, 0, 1, 2, ..., m
//      width of b[i] is M_VALUE, bit-order is little-endian
//   then:
//      px[POLY_X_WIDTH-1:0] = {b[k], ..., b[2], b[1], b[0]}
//
//////////////////////////////////////////////////////////////////////////////////

//
// Polynomial multiplication in GF(2)
// Operation: c(x) = a(x) * b(x)
//
module poly_mul #(
	parameter POLY_A_WIDTH = 4,		// Width of polynomial A
	parameter POLY_B_WIDTH = 4      // Width of polynomial A
)
(
    input wire [POLY_A_WIDTH-1:0]pa,					// Polynomial A
    input wire [POLY_B_WIDTH-1:0]pb,					// Polynomial B
    output wire [(POLY_A_WIDTH+POLY_B_WIDTH-1)-1:0]pc   // Polynomial C (product)
);

reg [(POLY_A_WIDTH+POLY_B_WIDTH-1)-1:0]and_table[POLY_A_WIDTH-1:0];
reg [(POLY_A_WIDTH+POLY_B_WIDTH-1)-1:0]xor_res;

integer i, j;

assign pc = xor_res;

always @(*) begin
	xor_res = {(POLY_A_WIDTH+POLY_B_WIDTH-1){1'b0}};
	
	for (i = 0; i < POLY_A_WIDTH; i = i + 1) begin
		and_table[i] = {((POLY_A_WIDTH+POLY_B_WIDTH)-1){1'b0}};
	end
	
	for (i = 0; i < POLY_A_WIDTH; i = i + 1) begin
		for (j = 0; j < POLY_B_WIDTH; j = j + 1) begin
			and_table[i][j+i] = pa[i] & pb[j];
		end
	end
	
	for (i = 0; i < POLY_A_WIDTH; i = i + 1) begin
		xor_res = xor_res ^ and_table[i];
	end
end

endmodule

//
// Polynomial division in GF(2)
// Operation: a(x) / b(x) -> a(x) = b(x) * q(x) + r(x)
//
module poly_div #(
	parameter POLY_A_WIDTH = 4,		// Width of polynomial A
	parameter POLY_B_WIDTH = 4		// Width of polynomial B
)
(
	input wire [POLY_A_WIDTH-1:0]pa,	// Polynomial A
	input wire [POLY_B_WIDTH-1:0]pb,	// Polynomial B
	output wire [POLY_A_WIDTH-1:0]pq,	// Polynomial Q (quotient) 
	output wire [POLY_B_WIDTH-1:0]pr	// Polynomial R (remainder) 
);

reg [(POLY_B_WIDTH-1)-1:0]r_state_tmp[POLY_A_WIDTH-1:0];
reg [POLY_B_WIDTH-1:0]c_state_tmp[POLY_A_WIDTH-1:0];
reg [(POLY_B_WIDTH-1)-1:0]r_state[POLY_A_WIDTH:0];
reg [POLY_A_WIDTH-1:0]quo;
reg [POLY_B_WIDTH-1:0]pb_deg;

integer i, j;

assign pr = r_state[POLY_A_WIDTH];
assign pq = quo;

always @(*) begin
	pb_deg = {POLY_B_WIDTH{1'b0}};
	
	for (i = 0; i < POLY_B_WIDTH; i = i + 1)begin
		if (pb[i] & 1'b1) begin
			pb_deg = i;
		end	 
	end
end

always @(*) begin
	r_state[0] = {(POLY_B_WIDTH-1){1'b0}};
	for (i = 0; i < POLY_A_WIDTH; i = i + 1) begin
		c_state_tmp[i][0] = pa[POLY_A_WIDTH-1-i];
		for (j = 1; j < POLY_B_WIDTH; j = j + 1) begin
			c_state_tmp[i][j] = r_state[i][j-1];
		end
		for (j = 0; j < POLY_B_WIDTH-1; j = j + 1) begin
			r_state_tmp[i][j] = c_state_tmp[i][j] ^ (pb[j] & c_state_tmp[i][pb_deg]);
		end
		quo[POLY_A_WIDTH-1-i] = c_state_tmp[i][pb_deg];
		r_state[i+1] = r_state_tmp[i];
	end
end

endmodule

//
// Polynomial modular multiplication in GF(2)
// Operation: c(x) = a(x)*b(x) (mod p(x))
//
module poly_mul_mod #(
	parameter POLY_A_WIDTH = 3,		// Width of Polynomial A
	parameter POLY_B_WIDTH = 2,		// Width of Polynomial B
	parameter POLY_P_WIDTH = 3		// Width of Polynomial P
)
(
	input wire [POLY_A_WIDTH-1:0]pa, 	// Polynomial A
	input wire [POLY_B_WIDTH-1:0]pb, 	// Polynomial B
	input wire [POLY_P_WIDTH-1:0]pp, 	// Polynomial P (must be irreducible over GF(2))
	output wire [POLY_P_WIDTH-1:0]pc	// Polynomial C (product)
);

wire [(POLY_A_WIDTH+POLY_B_WIDTH-1)-1:0]mul_res;
wire [POLY_P_WIDTH-1:0]rem_res;

assign pc = rem_res;
			
poly_mul #(
	.POLY_A_WIDTH(POLY_A_WIDTH),
	.POLY_B_WIDTH(POLY_B_WIDTH)
) poly_mul_inst (
	.pa(pa),
	.pb(pb),
	.pc(mul_res)
);
			
poly_div #(
	.POLY_A_WIDTH(POLY_A_WIDTH+POLY_B_WIDTH-1),
	.POLY_B_WIDTH(POLY_P_WIDTH)
) poly_div_inst (
	.pa(mul_res),
	.pb(pp),
	.pq(),
	.pr(rem_res)
);

endmodule

//
// Polynomial multiplication in GF(2^m)
// Operation: c(x) = a(x) * b(x)
//
module poly_mul_2m #(
	parameter POLY_A_WIDTH = 3,		// Width of Polynomial A
	parameter POLY_B_WIDTH = 2,		// Width of Polynomial B
	parameter M_VALUE = 2			// Degree of 2 (m)
)
(
	input wire [POLY_A_WIDTH*M_VALUE-1:0]pa,					// Polynomial A
	input wire [POLY_B_WIDTH*M_VALUE-1:0]pb,					// Polynomial B
	input wire [M_VALUE:0]pp,									// Polynomial P (must be irreducible over GF(2))
	output wire [(POLY_A_WIDTH+POLY_B_WIDTH-1)*M_VALUE-1:0]pc	// Polynomial C (product)
);

reg [(POLY_A_WIDTH+POLY_B_WIDTH-1)*M_VALUE-1:0]xor_table[POLY_A_WIDTH-1:0];
reg [(POLY_A_WIDTH+POLY_B_WIDTH-1)*M_VALUE-1:0]xor_res;
wire [POLY_B_WIDTH*M_VALUE-1:0]mul_mod_table[POLY_A_WIDTH-1:0];

integer i;
genvar g, k;

assign pc = xor_res;

always @(*) begin
	xor_res = {((POLY_A_WIDTH+POLY_B_WIDTH-1)*M_VALUE){1'b0}};
	
	for (i = 0; i < POLY_A_WIDTH; i = i + 1) begin
		xor_table[i] = {((POLY_A_WIDTH+POLY_B_WIDTH-1)*M_VALUE){1'b0}};
	end
	
	for (i = 0; i < POLY_A_WIDTH; i = i + 1) begin
		xor_table[i][(POLY_B_WIDTH+i)*M_VALUE-1-:POLY_B_WIDTH*M_VALUE] = mul_mod_table[i];
	end
	
	for (i = 0; i < POLY_A_WIDTH; i = i + 1) begin
		xor_res = xor_res ^ xor_table[i];
	end
end

generate
	for (g = 0; g < POLY_A_WIDTH; g = g + 1) begin
		for (k = 0; k < POLY_B_WIDTH; k = k + 1) begin
			wire [M_VALUE:0]mul_mod_res;
			
			assign mul_mod_table[g][(k+1)*M_VALUE-1-:M_VALUE] = mul_mod_res[M_VALUE-1:0];
			
			poly_mul_mod #(
				.POLY_A_WIDTH(M_VALUE),
				.POLY_B_WIDTH(M_VALUE),
				.POLY_P_WIDTH(M_VALUE+1)
			) poly_mul_mod_inst (
				.pa(pa[(g+1)*M_VALUE-1-:M_VALUE]),
				.pb(pb[(k+1)*M_VALUE-1-:M_VALUE]),
				.pp(pp),
				.pc(mul_mod_res)
			);
			
		end
	end
endgenerate

endmodule

//
// Polynomial division in GF(2^m)
// Operation: a(x) / b(x) -> a(x) = b(x) * q(x) + r(x)
// b(x) must be type of b[0]+b[1]*x+b[2]*x^2+...+b[k-1]*x^(k-1)+1*x^k, where k <= (POLY_B_WIDTH-1) 
//   and b[i] is element of GF(2^m)
//
module poly_div_2m #(
	parameter POLY_A_WIDTH = 3,		// Width of Polynomial A
	parameter POLY_B_WIDTH = 2,		// Width of Polynomial B
	parameter M_VALUE = 2			// Degree of 2 (m)
)
(
	input wire [POLY_A_WIDTH*M_VALUE-1:0]pa,	// Polynomial A
	input wire [POLY_B_WIDTH*M_VALUE-1:0]pb,	// Polynomial B
	input wire [M_VALUE:0]pp,					// Polynomial P (must be irreducible over GF(2))
	output wire [POLY_A_WIDTH*M_VALUE-1:0]pq,	// Polynomial Q (quotient) 
	output wire [POLY_B_WIDTH*M_VALUE-1:0]pr	// Polynomial R (remainder) 
);

reg [POLY_B_WIDTH*M_VALUE-1:0]c_state[POLY_A_WIDTH-1:0];
reg [(POLY_B_WIDTH-1)*M_VALUE-1:0]r_state[POLY_A_WIDTH:0];
reg [(POLY_B_WIDTH-1)*M_VALUE-1:0]r_state_tmp[POLY_A_WIDTH-1:0];
reg [POLY_A_WIDTH*M_VALUE-1:0]quo;
reg [POLY_B_WIDTH*M_VALUE-1:0]pb_deg;

wire [(POLY_B_WIDTH-1)*M_VALUE-1:0]mul_2m_table[POLY_A_WIDTH-1:0];

integer i, j;
genvar g, k;

assign pr = r_state[POLY_A_WIDTH];
assign pq = quo;

always @(*) begin
	pb_deg = {(POLY_B_WIDTH*M_VALUE){1'b0}};
	
	for (i = 0; i < POLY_B_WIDTH; i = i + 1)begin
		if (pb[(i+1)*M_VALUE-1-:M_VALUE] != {(M_VALUE){1'b0}}) begin
			pb_deg = i;
		end	 
	end
end

always @(*) begin
	r_state[0] = {((POLY_B_WIDTH-1)*M_VALUE){1'b0}};
	
	for (i = 0; i < POLY_A_WIDTH; i = i + 1) begin
		c_state[i][M_VALUE-1:0] = pa[(POLY_A_WIDTH-i)*M_VALUE-1-:M_VALUE];
		for (j = 1; j < POLY_B_WIDTH; j = j + 1) begin
			c_state[i][(j+1)*M_VALUE-1-:M_VALUE] = r_state[i][j*M_VALUE-1-:M_VALUE];
		end
		for (j = 0; j < POLY_B_WIDTH-1; j = j + 1) begin
			r_state_tmp[i][(j+1)*M_VALUE-1-:M_VALUE] = 
			c_state[i][(j+1)*M_VALUE-1-:M_VALUE] ^ mul_2m_table[i][(j+1)*M_VALUE-1-:M_VALUE];
		end
		quo[(POLY_A_WIDTH-i)*M_VALUE-1-:M_VALUE] = c_state[i][(pb_deg+1)*M_VALUE-1-:M_VALUE];
		r_state[i+1] = r_state_tmp[i];
	end
end

generate
	for (g = 0; g < POLY_A_WIDTH; g = g + 1) begin
		for (k = 0; k < POLY_B_WIDTH - 1; k = k + 1) begin

			wire [M_VALUE:0]mul_mod_res;

			assign mul_2m_table[g][(k+1)*M_VALUE-1-:M_VALUE] = mul_mod_res[M_VALUE-1:0];
			
			poly_mul_mod #(
				.POLY_A_WIDTH(M_VALUE),
				.POLY_B_WIDTH(M_VALUE),
				.POLY_P_WIDTH(M_VALUE+1)
			) poly_mul_mod_inst (
				.pa(c_state[g][(pb_deg+1)*M_VALUE-1-:M_VALUE]),
				.pb(pb[(k+1)*M_VALUE-1-:M_VALUE]),
				.pp(pp),
				.pc(mul_mod_res)
			);
		
		end
	end
endgenerate

endmodule

//
// Greatest Common Divisor of two polynomials in GF(2)
// Operation: c(x) = GCD(a(x), b(x))
//
module poly_gcd #(
	parameter POLY_WIDTH = 4 	// Width of Polynomials
)
(
	input wire [POLY_WIDTH-1:0]pa,		// Polynomial A
	input wire [POLY_WIDTH-1:0]pb,		// Polynomial B
	output wire [POLY_WIDTH-1:0]pc		// Polynomial C (product)
);

wire [POLY_WIDTH-1:0]pa_vec[POLY_WIDTH:0];
wire [POLY_WIDTH-1:0]pb_vec[POLY_WIDTH:0];
wire [POLY_WIDTH-1:0]pr_vec[POLY_WIDTH-1:0];
reg [POLY_WIDTH-1:0]gcd;

integer i;
genvar g;

assign pa_vec[0] = (pa > pb) ? pa : pb;
assign pb_vec[0] = (pa > pb) ? pb : pa;

assign pc = gcd;

always @(*) begin
	gcd = {{(POLY_WIDTH-1){1'b0}}, {1'b1}};
	
	for (i = POLY_WIDTH - 1; i >= 0; i = i - 1) begin
		if (pr_vec[i] == {POLY_WIDTH{1'b0}}) begin
			gcd = pb_vec[i];
		end
	end
end

generate
	for (g = 0; g < POLY_WIDTH; g = g + 1) begin
		
		assign pa_vec[g+1] = pb_vec[g];
		assign pb_vec[g+1] = pr_vec[g];
	
		poly_div #(
			.POLY_A_WIDTH(POLY_WIDTH),
			.POLY_B_WIDTH(POLY_WIDTH)
		) poly_width_inst (
			.pa(pa_vec[g]),
			.pb(pb_vec[g]),
			.pq(),
			.pr(pr_vec[g])
		);
	end
endgenerate

endmodule

//
// Least Common Multiple of two polynomials in GF(2)
// Operation: c(x) = LCM(a(x), b(x))
//
module poly_lcm #(
	parameter POLY_WIDTH = 4	// Width of Polynomials
)
( 
	input wire [POLY_WIDTH-1:0]pa,			// Polynomial A
	input wire [POLY_WIDTH-1:0]pb,			// Polynomial B
	output wire [(2*POLY_WIDTH-1)-1:0]pc	// Polynomial C (product)
);

wire [POLY_WIDTH-1:0]gcd;
wire [(2*POLY_WIDTH-1)-1:0]mul;

poly_gcd #(
	.POLY_WIDTH(POLY_WIDTH)
) poly_gcd_inst (
	.pa(pa),
	.pb(pb),
	.pc(gcd)
);

poly_mul #(
	.POLY_A_WIDTH(POLY_WIDTH),
	.POLY_B_WIDTH(POLY_WIDTH)
) poly_mul_inst (
    .pa(pa),
    .pb(pb),
    .pc(mul)
);

poly_div #(
	.POLY_A_WIDTH(2*POLY_WIDTH-1),
	.POLY_B_WIDTH(POLY_WIDTH)
) poly_div_inst (
	.pa(mul),
	.pb(gcd),
	.pq(pc),
	.pr()
);

endmodule